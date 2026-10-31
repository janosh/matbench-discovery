"""DPA4/DPA4C export routing, device guard and kernel cache separation."""

# ruff: noqa: SLF001

import sys
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from matbench_discovery import calculators


@pytest.fixture
def fake_deepmd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Stub torch and DeePMD so export policy runs without a GPU."""
    gpu = SimpleNamespace(type="cuda", index=0)
    cuda = SimpleNamespace(get_device_name=lambda _device: "NVIDIA A800-SXM4-80GB")
    torch = SimpleNamespace(
        device=lambda name: (
            gpu if name == "cuda" else SimpleNamespace(type=name, index=None)
        ),
        cuda=cuda,
    )
    calc_cls, exporter = Mock(), Mock()
    monkeypatch.setitem(sys.modules, "torch", torch)
    monkeypatch.setitem(sys.modules, "deepmd.calculator", SimpleNamespace(DP=calc_cls))
    monkeypatch.setitem(
        sys.modules, "deepmd.pt_expt.utils.env", SimpleNamespace(DEVICE=gpu)
    )
    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", str(tmp_path))
    monkeypatch.setattr(calculators, "_run_to_atomic_output", exporter)
    monkeypatch.delenv("DP_TRITON_INFER", raising=False)
    monkeypatch.delenv("DP_CUDA_INFER", raising=False)
    checkpoint = tmp_path / "model.pt"
    checkpoint.write_bytes(b"weights")
    return SimpleNamespace(
        checkpoint=str(checkpoint), calc_cls=calc_cls, exporter=exporter, cuda=cuda
    )


@pytest.mark.parametrize(
    ("factory", "backend", "options"),
    [
        (calculators._DPA4_GRAPH, "--pt", {"nlist_backend": "native"}),
        (calculators._DPA4C_GRAPH, "--pt-expt", {"neighbor_graph_method": "cell"}),
    ],
)
def test_deepmd_graph_export_and_load(
    fake_deepmd: SimpleNamespace,
    factory: Callable[..., object],
    backend: str,
    options: dict[str, str],
) -> None:
    """Each architecture exports with its CLI and loads with its neighbor builder."""
    factory("cuda", fake_deepmd.checkpoint)
    (command, frozen), kwargs = fake_deepmd.exporter.call_args
    checkpoint = fake_deepmd.checkpoint
    assert command == ["dp", backend, "freeze", "-c", checkpoint, "-o", "{output}"]
    assert frozen.endswith(f"-{backend[2:]}-NVIDIA-A800-SXM4-80GB-2-1.pt2")
    assert kwargs["env_overrides"] == {"DP_TRITON_INFER": "2", "DP_CUDA_INFER": "1"}
    fake_deepmd.calc_cls.assert_called_once_with(frozen, **options)

    # an existing .pt2 skips export
    fake_deepmd.exporter.reset_mock()
    factory("cuda", "validated.pt2")
    fake_deepmd.exporter.assert_not_called()
    fake_deepmd.calc_cls.assert_called_with("validated.pt2", **options)


def test_deepmd_graph_cache_key_tracks_gpu_and_kernels(
    fake_deepmd: SimpleNamespace, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A different GPU model or kernel policy never reuses a previous export."""
    paths = []
    for gpu_name, triton in [("A800", "2"), ("A100", "2"), ("A100", "0")]:
        monkeypatch.setattr(
            fake_deepmd.cuda, "get_device_name", lambda _, n=gpu_name: n
        )
        monkeypatch.setenv("DP_TRITON_INFER", triton)
        calculators._DPA4_GRAPH("cuda", fake_deepmd.checkpoint)
        paths.append(fake_deepmd.exporter.call_args.args[1])
    assert len(set(paths)) == 3
    assert paths[2].endswith("-A100-0-1.pt2")


def test_deepmd_graph_rejects_inconsistent_device(fake_deepmd: SimpleNamespace) -> None:
    """A CPU request must not silently run on the CUDA device DeePMD imported with."""
    with pytest.raises(ValueError, match="DeePMD initialized on"):
        calculators._DPA4_GRAPH("cpu", fake_deepmd.checkpoint)
    fake_deepmd.exporter.assert_not_called()
    fake_deepmd.calc_cls.assert_not_called()
