"""DPA4/DPA4C export routing, device guard and kernel cache separation."""

# ruff: noqa: SLF001

import sys
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from matbench_discovery import calculators


def fake_device(name: str) -> SimpleNamespace:
    """Mimic torch.device: 'cuda' has index None, 'cuda:1' index 1."""
    dev_type, _, index = name.partition(":")
    return SimpleNamespace(type=dev_type, index=int(index) if index else None)


@pytest.fixture
def fake_deepmd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Stub torch and DeePMD (imported on cuda:0) so export policy runs on CPU."""
    cuda = SimpleNamespace(get_device_name=Mock(return_value="NVIDIA A800 (SXM4)"))
    env = SimpleNamespace(DEVICE=fake_device("cuda:0"))
    calc_cls, exporter = Mock(), Mock()
    torch = SimpleNamespace(device=fake_device, cuda=cuda)
    monkeypatch.setitem(sys.modules, "torch", torch)
    monkeypatch.setitem(sys.modules, "deepmd.calculator", SimpleNamespace(DP=calc_cls))
    monkeypatch.setitem(sys.modules, "deepmd.pt_expt.utils.env", env)
    monkeypatch.setattr(calculators, "CHECKPOINT_DIR", str(tmp_path))
    monkeypatch.setattr(calculators, "_run_to_atomic_output", exporter)
    for key in ("DP_TRITON_INFER", "DP_CUDA_INFER", "DP_AMP_INFER"):
        monkeypatch.delenv(key, raising=False)
    checkpoint = tmp_path / "model.pt"
    checkpoint.write_bytes(b"weights")
    return SimpleNamespace(
        checkpoint=str(checkpoint),
        calc_cls=calc_cls,
        exporter=exporter,
        cuda=cuda,
        env=env,
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
    assert frozen.endswith(f"-{backend[2:]}-NVIDIA-A800-SXM4-2-1-0.pt2")
    assert kwargs["env_overrides"] == {
        "DP_TRITON_INFER": "2",
        "DP_CUDA_INFER": "1",
        "DP_AMP_INFER": "0",
    }
    assert kwargs["timeout_sec"] == 3 * 60 * 60
    fake_deepmd.calc_cls.assert_called_once_with(frozen, **options)

    # an existing .pt2 skips export
    fake_deepmd.exporter.reset_mock()
    factory("cuda", "validated.pt2")
    fake_deepmd.exporter.assert_not_called()
    fake_deepmd.calc_cls.assert_called_with("validated.pt2", **options)


def test_deepmd_graph_cache_key_tracks_inputs(
    fake_deepmd: SimpleNamespace, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Weights, GPU model and kernel/precision policy each get their own export."""
    other_checkpoint = tmp_path / "other.pt"
    other_checkpoint.write_bytes(b"other weights")
    variants = [
        (fake_deepmd.checkpoint, "A800", {}),
        (str(other_checkpoint), "A800", {}),
        (fake_deepmd.checkpoint, "A100", {}),
        (fake_deepmd.checkpoint, "A100", {"DP_TRITON_INFER": "0"}),
        (fake_deepmd.checkpoint, "A100", {"DP_AMP_INFER": "1"}),
    ]
    paths = []
    for checkpoint, gpu_name, env in variants:
        fake_deepmd.cuda.get_device_name.return_value = gpu_name
        for key, value in env.items():
            monkeypatch.setenv(key, value)
        calculators._DPA4_GRAPH("cuda", checkpoint)
        paths.append(fake_deepmd.exporter.call_args.args[1])
        for key in env:
            monkeypatch.delenv(key)
    assert len(set(paths)) == len(variants)


@pytest.mark.parametrize(
    ("requested", "imported", "suffix"),
    [
        ("cpu", "cpu", "-cpu-0-0-0.pt2"),
        ("cuda:1", "cuda:0", None),
        ("cpu", "cuda:0", None),
    ],
)
def test_deepmd_graph_device_guard(
    fake_deepmd: SimpleNamespace, requested: str, imported: str, suffix: str | None
) -> None:
    """Matching devices export (CPU without GPU kernels); mismatches never run."""
    fake_deepmd.env.DEVICE = fake_device(imported)
    if suffix is None:
        with pytest.raises(ValueError, match="DeePMD initialized on"):
            calculators._DPA4_GRAPH(requested, fake_deepmd.checkpoint)
        fake_deepmd.exporter.assert_not_called()
        fake_deepmd.calc_cls.assert_not_called()
        return
    calculators._DPA4_GRAPH(requested, fake_deepmd.checkpoint)
    assert fake_deepmd.exporter.call_args.args[1].endswith(suffix)
    if requested == "cpu":
        fake_deepmd.cuda.get_device_name.assert_not_called()
