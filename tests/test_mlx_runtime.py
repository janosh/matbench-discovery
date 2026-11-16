"""Mixed-backend provenance must report actual MLX work instead of CPU weights."""

import sys
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from matbench_discovery import hpc


@pytest.fixture
def mlx(monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    class Array:
        def __init__(self, dtype: str = "float32") -> None:
            self.dtype = dtype

    gpu = object()
    module = SimpleNamespace(
        array=Array,
        gpu=gpu,
        default_device=lambda: gpu,
        device_info=lambda: {"device_name": "Apple Test GPU"},
        get_peak_memory=lambda: 750_000_000,
        reset_peak_memory=Mock(),
    )
    monkeypatch.setitem(sys.modules, "mlx.core", module)
    return module


def test_mlx_opt_in_precedes_torch_geometry_master_weights(
    mlx: SimpleNamespace, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        hpc, "_torch_runtime", lambda _: pytest.fail("wrong CPU attribution")
    )
    calc = SimpleNamespace(mbd_backend="mlx", mbd_mlx_arrays=[mlx.array()])
    assert hpc.detect_runtime(calc, requested_device="mps") == {
        "backend": "mlx",
        "device": "mps",
        "dtype": "float32",
        "hardware": "Apple Test GPU (Metal)",
    }


def test_installed_mlx_does_not_relabel_other_models(
    mlx: SimpleNamespace, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert sys.modules["mlx.core"] is mlx
    expected = {"backend": "torch", "device": "cpu", "hardware": "CPU"}
    monkeypatch.setattr(hpc, "_torch_runtime", lambda _: expected)
    assert hpc.detect_runtime(object(), requested_device="cpu") == expected


@pytest.mark.parametrize("problem", ["not_loaded", "cpu", "no_weights", "mixed_dtype"])
def test_mlx_provenance_fails_closed(
    mlx: SimpleNamespace, monkeypatch: pytest.MonkeyPatch, problem: str
) -> None:
    calc = SimpleNamespace(mbd_backend="mlx", mbd_mlx_arrays=[mlx.array()])
    if problem == "not_loaded":
        monkeypatch.delitem(sys.modules, "mlx.core")
    if problem == "cpu":
        mlx.default_device = lambda: "cpu"
    if problem == "no_weights":
        calc.mbd_mlx_arrays = []
    if problem == "mixed_dtype":
        calc.mbd_mlx_arrays.append(mlx.array("float16"))
    with pytest.raises(ValueError, match="MLX"):
        hpc.detect_runtime(calc, requested_device="mps")


def test_mlx_gpu_memory_is_measured_and_reset(mlx: SimpleNamespace) -> None:
    hpc.reset_gpu_peak_memory({"backend": "mlx"})
    mlx.reset_peak_memory.assert_called_once_with()
    result = hpc.memory_peaks({"backend": "mlx"}, succeeded=True)
    assert result["max_gpu_mem_gb"] == 0.75


def test_successful_mlx_run_requires_measured_gpu_allocation(
    mlx: SimpleNamespace,
) -> None:
    mlx.get_peak_memory = lambda: 0
    with pytest.raises(ValueError, match="MLX"):
        hpc.memory_peaks({"backend": "mlx"}, succeeded=True)


@pytest.mark.parametrize("device", ["cpu", "cuda", "mps"])
def test_shared_sharded_runner_accepts_schema_supported_device(device: str) -> None:
    import argparse

    from matbench_discovery.runner_cli import add_shard_args

    parser = argparse.ArgumentParser()
    add_shard_args(parser, merge_shards_help="merge", write_yaml_help="write")
    assert parser.parse_args(["--device", device]).device == device
