"""Tests for high-performance computing utilities."""

import os
import sys
from collections.abc import Sequence
from unittest.mock import mock_open, patch

import pandas as pd
import pytest

from matbench_discovery import hpc


@pytest.mark.parametrize(
    "n_chunks,task_id,expected_index",
    [(3, 1, [0, 1, 2, 3]), (3, 2, [4, 5, 6]), (3, 3, [7, 8, 9]), (1, 1, range(10))],
)
def test_df_slurm_chunk(
    n_chunks: int, task_id: int, expected_index: Sequence[int]
) -> None:
    """df_slurm_chunk returns the right rows for a 1-based slurm array task."""
    df_in = pd.DataFrame({"col": range(10)})
    df_chunk = hpc.df_slurm_chunk(df_in, n_chunks, task_id)
    assert isinstance(df_chunk, pd.DataFrame)
    assert list(df_chunk.index) == list(expected_index)


@pytest.mark.parametrize("n_chunks,task_id", [(0, 1), (3, 0), (3, -1), (3, 4)])
def test_df_slurm_chunk_rejects_invalid_task_ids(n_chunks: int, task_id: int) -> None:
    """df_slurm_chunk rejects invalid chunk counts and task IDs."""
    df_in = pd.DataFrame({"col": range(10)})
    with pytest.raises(ValueError, match=rf"{n_chunks=}.*{task_id=}"):
        hpc.df_slurm_chunk(df_in, n_chunks, task_id)


@pytest.mark.parametrize("py_file_path", [None, "path/to/file.py"])
@pytest.mark.parametrize("partition", [None, "fake-partition"])
@pytest.mark.parametrize("time", [None, "0:0:1"])
@pytest.mark.parametrize("account", [None, "fake-account"])
@pytest.mark.parametrize("pre_cmd", [None, "module load pytorch;", "ENV_VAR=42"])
@pytest.mark.parametrize("submit_as_temp_file", [True, False])
def test_slurm_submit(
    capsys: pytest.CaptureFixture[str],
    py_file_path: str | None,
    partition: str | None,
    time: str | None,
    account: str | None,
    pre_cmd: str | None,
    submit_as_temp_file: bool,
) -> None:
    job_name = "test_job"
    out_dir = "tmp"

    kwargs = dict(
        job_name=job_name,
        out_dir=out_dir,
        time=time,
        partition=partition,
        account=account,
        py_file_path=py_file_path,
        slurm_flags="--foo",
        pre_cmd=pre_cmd,
        submit_as_temp_file=submit_as_temp_file,
    )

    hpc.slurm_submit(**kwargs)  # ty: ignore[invalid-argument-type]

    stdout, stderr = capsys.readouterr()
    # check slurm_submit() did nothing in normal mode
    assert stdout == stderr == ""

    with patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True):
        slurm_vars = hpc.slurm_submit(**kwargs)  # ty: ignore[invalid-argument-type]
    expected_slurm_vars = dict(slurm_job_id="1234", slurm_flags="--foo")
    if time is not None:
        expected_slurm_vars["slurm_timelimit"] = time
    if pre_cmd and not pre_cmd.strip().endswith(";"):
        pre_cmd += ";"
    if pre_cmd is not None and pre_cmd != "":
        expected_slurm_vars["pre_cmd"] = pre_cmd
    assert slurm_vars == expected_slurm_vars

    # check slurm_submit() prints cmd and calls subprocess.run() in submit mode
    with (
        pytest.raises(SystemExit),
        patch("sys.argv", ["slurm-submit"]),
        patch("matbench_discovery.hpc.subprocess.run") as mock_subprocess_run,
        patch(
            "matbench_discovery.hpc.tempfile.mkdtemp",
            return_value="/tmp/slurm_job_123",
        ),
        patch("matbench_discovery.hpc.shutil.copy2") as mock_copy2,
        patch("builtins.open", mock_open()),
        patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True),
    ):
        hpc.slurm_submit(**kwargs)  # ty: ignore[invalid-argument-type]

    assert mock_subprocess_run.call_count == 1

    expected_py_file_path = py_file_path or __file__
    if submit_as_temp_file:
        expected_py_file_path = (
            f"/tmp/slurm_job_123/{os.path.basename(expected_py_file_path)}"
        )
        assert mock_copy2.called
        assert mock_copy2.call_args[0][0] == (py_file_path or __file__)
        assert mock_copy2.call_args[0][1] == expected_py_file_path
    else:
        assert not mock_copy2.called

    sbatch_cmd = (
        f"sbatch --job-name {job_name} --output {out_dir}/slurm-%A.log --foo "
        f"--wrap {pre_cmd + ' ' if pre_cmd else ''}python {expected_py_file_path}"
    ).replace(" --", "\n  --")
    for flag in (f"{time=!s}", f"{account=!s}", f"{partition=!s}"):
        key, val = flag.split("=")
        if val != "None":
            sbatch_cmd += f"\n  --{key} {val}"

    stdout, stderr = capsys.readouterr()
    assert sbatch_cmd in stdout
    assert stderr == ""


def test_get_calling_file_path() -> None:
    assert hpc._get_calling_file_path(frame=1) == __file__  # noqa: SLF001

    def wrapper(frame: int) -> str:
        return hpc._get_calling_file_path(frame)  # noqa: SLF001

    assert wrapper(frame=2) == __file__


@pytest.mark.parametrize(
    ("shard_metadatas", "expected"),
    [
        (
            [
                # unrelated keys like excluded_formula_reasons are ignored by merge
                {
                    "hardware": "H200",
                    "run_time_sec": 100.0,
                    "excluded_formula_reasons": {},
                },
                {"hardware": "H200", "run_time_sec": 50.5},
            ],
            {"hardware": "H200", "run_time_sec": 150.5},
        ),
        ([{}, {}], {}),  # no provenance -> empty (leaves existing YAML untouched)
        ([{"run_time_sec": 12.0}, {}], {"run_time_sec": 12.0}),  # partial
        (
            # memory peaks max over shards (each shard is its own process)
            [
                {"run_time_sec": 10.0, "max_rss_gb": 4.2, "max_gpu_mem_gb": 11.5},
                {"run_time_sec": 20.0, "max_rss_gb": 6.1, "max_gpu_mem_gb": 9.0},
                {"run_time_sec": 5.0},  # shard without memory info still sums time
            ],
            {"run_time_sec": 35.0, "max_rss_gb": 6.1, "max_gpu_mem_gb": 11.5},
        ),
    ],
    ids=["full", "missing", "partial", "memory_peaks"],
)
def test_merge_run_metadata(
    shard_metadatas: list[dict[str, object]], expected: dict[str, str | float]
) -> None:
    """merge_run_metadata sums run_time_sec, maxes memory peaks + shares hardware."""
    assert hpc.merge_run_metadata(shard_metadatas) == expected


INCOMPLETE_HARDWARE_SEGMENTS = (
    {"hardware": "H200", "versions": {}},
    {"hardware": "", "versions": {}},
)


def test_merge_audit_metadata_drops_partial_cost_fields() -> None:
    """Audit aggregation publishes each cost field only with complete coverage."""
    merged = hpc.merge_audit_metadata(
        [
            {"hardware": "H200", "run_time_sec": 10.0, "max_rss_gb": 4.0},
            {"hardware": "H200", "run_time_sec": 20.0},
        ]
    )
    assert merged == {"hardware": "H200", "run_time_sec": 30.0}
    base_segment = {"versions": {"numpy": "1"}}
    for other_segment in ({}, {"versions": {"numpy": "2"}}):
        merged = hpc.merge_audit_metadata([base_segment, other_segment])
        assert "versions" not in merged


def test_merge_audit_metadata_drops_mismatched_hardware() -> None:
    """Mismatched hardware provenance is silently omitted from merged metadata."""
    assert "hardware" not in hpc.merge_audit_metadata(INCOMPLETE_HARDWARE_SEGMENTS)


def test_merge_audit_metadata_strictly_rejects_missing_hardware() -> None:
    """Strict audit merging rejects a segment with missing hardware provenance."""
    with pytest.raises(ValueError, match="missing hardware"):
        hpc.merge_audit_metadata(INCOMPLETE_HARDWARE_SEGMENTS, strict=True)


def test_peak_memory_gb() -> None:
    """peak_memory_gb reports a positive Unix RSS peak; GPU reset is a safe no-op."""
    mem = hpc.peak_memory_gb()
    if sys.platform != "win32":  # Windows lacks the resource module
        assert mem["max_rss_gb"] > 0
    assert set(mem) <= {"max_rss_gb", "max_gpu_mem_gb"}
    hpc.reset_gpu_peak_memory()  # must not raise without torch/CUDA
