"""DPA4/DPA4C export routing, device guards, and kernel cache separation."""

# ruff: noqa: SLF001

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from ase.calculators.calculator import Calculator

from matbench_discovery import calculators


@pytest.fixture
def fake_deepmd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Isolate export policy tests from actual GPU compilation."""
    device = SimpleNamespace(type="cuda", index=0)
    torch = SimpleNamespace(
        device=lambda value: (
            device if value == "cuda" else SimpleNamespace(type=value, index=None)
        ),
        cuda=SimpleNamespace(
            get_device_capability=lambda _device: (8, 0),
            get_device_name=lambda _device: "NVIDIA A800-SXM4-80GB",
        ),
    )
    model = SimpleNamespace(get_type_map=lambda: ["H", "O"], get_ntypes=lambda: 2)
    loader = Mock(return_value=model)
    exporter = Mock()
    monkeypatch.setitem(sys.modules, "torch", torch)
    monkeypatch.setitem(
        sys.modules, "deepmd.calculator", SimpleNamespace(DP=Calculator)
    )
    monkeypatch.setitem(sys.modules, "deepmd.infer", SimpleNamespace(DeepPot=loader))
    monkeypatch.setitem(
        sys.modules, "deepmd.pt_expt.utils.env", SimpleNamespace(DEVICE=device)
    )
    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", str(tmp_path / "cache"))
    monkeypatch.setattr(calculators, "_run_to_atomic_output", exporter)
    monkeypatch.delenv("DP_TRITON_INFER", raising=False)
    monkeypatch.delenv("DP_CUDA_INFER", raising=False)
    checkpoint = tmp_path / "model.pt"
    checkpoint.write_bytes(b"distinct model weights")
    return SimpleNamespace(checkpoint=str(checkpoint), loader=loader, exporter=exporter)


@pytest.mark.parametrize("backend", ["--pt", "--pt-expt"])
def test_deepmd_graph_routes_backend_and_neighbor_interface(
    fake_deepmd: SimpleNamespace, backend: str
) -> None:
    """Each architecture gets its proper export CLI and unchanged DP calculation."""
    calc = calculators._deepmd_graph(backend)("cuda", fake_deepmd.checkpoint)
    command, output = fake_deepmd.exporter.call_args.args
    assert command == [
        "dp",
        backend,
        "freeze",
        "-c",
        fake_deepmd.checkpoint,
        "-o",
        "{output}",
    ]
    assert f"-{backend[2:]}-sm80-" in output
    assert output.endswith("-2-1.pt2")
    expected_options = (
        {"nlist_backend": "native"}
        if backend == "--pt"
        else {"neighbor_graph_method": "cell"}
    )
    fake_deepmd.loader.assert_called_once_with(output, **expected_options)
    assert calc.type_dict == {"H": 0, "O": 1}
    assert fake_deepmd.exporter.call_args.kwargs["tool_packages"] == (
        "deepmd-kit",
        "torch",
        "triton",
    )


def test_deepmd_graph_cache_separates_kernel_policy(
    fake_deepmd: SimpleNamespace, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An export with different baked kernels cannot reuse the previous archive."""
    factory = calculators._deepmd_graph("--pt")
    factory("cuda", fake_deepmd.checkpoint)
    previous = fake_deepmd.exporter.call_args.args[1]
    monkeypatch.setenv("DP_TRITON_INFER", "0")
    factory("cuda", fake_deepmd.checkpoint)
    current = fake_deepmd.exporter.call_args.args[1]
    assert previous != current
    assert "-sm80-" in current
    assert current.endswith("-0-1.pt2")


def test_deepmd_graph_cache_separates_gpu_models(
    fake_deepmd: SimpleNamespace, monkeypatch: pytest.MonkeyPatch
) -> None:
    """GPU models sharing a compute capability still have distinct tuning tables."""
    factory = calculators._deepmd_graph("--pt")
    factory("cuda", fake_deepmd.checkpoint)
    previous = fake_deepmd.exporter.call_args.args[1]
    monkeypatch.setattr(
        sys.modules["torch"].cuda,
        "get_device_name",
        lambda _device: "NVIDIA A100-SXM4-80GB",
    )
    factory("cuda", fake_deepmd.checkpoint)
    assert previous != fake_deepmd.exporter.call_args.args[1]


def test_deepmd_graph_accepts_frozen_override(fake_deepmd: SimpleNamespace) -> None:
    """Already exported local artifacts avoid compilation and downloading."""
    calculators._deepmd_graph("--pt-expt")("cuda", "validated.pt2")
    fake_deepmd.exporter.assert_not_called()
    fake_deepmd.loader.assert_called_once_with(
        "validated.pt2", neighbor_graph_method="cell"
    )


def test_deepmd_graph_rejects_inconsistent_device(fake_deepmd: SimpleNamespace) -> None:
    """A CPU request must not silently execute on the imported CUDA backend."""
    with pytest.raises(ValueError, match="DeePMD initialized on"):
        calculators._deepmd_graph("--pt")("cpu", fake_deepmd.checkpoint)
    fake_deepmd.exporter.assert_not_called()
    fake_deepmd.loader.assert_not_called()
