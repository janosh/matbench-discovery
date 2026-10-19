"""Tests for high-performance computing utilities."""

import os
import sys
from collections.abc import Sequence
from pathlib import Path
from unittest.mock import patch

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


@pytest.mark.parametrize(
    (
        "py_file_path",
        "time",
        "account",
        "partition",
        "array",
        "pre_cmd",
        "submit_as_temp_file",
    ),
    [
        (None, None, None, None, None, "", True),
        (
            "path/to/file.py",
            "0:0:1",
            "fake-account",
            "fake-partition",
            None,
            "module load pytorch;",
            False,
        ),
        (None, "0:0:1", None, None, "0-9%2", "ENV_VAR=42", True),
        ("path/to/file.py", None, "fake-account", None, "1-10", "", False),
        (None, None, None, "fake-partition", None, "ENV_VAR=42", False),
    ],
)
def test_slurm_submit(
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
    py_file_path: str | None,
    time: str | None,
    account: str | None,
    partition: str | None,
    array: str | None,
    pre_cmd: str,
    submit_as_temp_file: bool,
) -> None:
    """slurm_submit is a no-op outside submit mode and builds the sbatch command."""
    job_name, out_dir = "test_job", f"{tmp_path}/slurm-logs"
    kwargs = dict(
        job_name=job_name,
        out_dir=out_dir,
        time=time,
        partition=partition,
        account=account,
        py_file_path=py_file_path,
        slurm_flags="--foo",
        array=array,
        pre_cmd=pre_cmd,
        submit_as_temp_file=submit_as_temp_file,
    )

    hpc.slurm_submit(**kwargs)  # ty: ignore[invalid-argument-type]
    assert capsys.readouterr() == ("", "")  # does nothing in normal mode
    assert os.path.isdir(out_dir)  # sbatch fails if the log dir is missing

    with patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True):
        slurm_vars = hpc.slurm_submit(**kwargs)  # ty: ignore[invalid-argument-type]
    expected_pre_cmd = (
        pre_cmd if pre_cmd.endswith(";") or not pre_cmd else f"{pre_cmd};"
    )
    assert slurm_vars == {
        "slurm_job_id": "1234",
        "slurm_flags": "--foo",
        **({"slurm_timelimit": time} if time else {}),
        **({"pre_cmd": expected_pre_cmd} if pre_cmd else {}),
    }

    # submit mode prints the cmd and calls subprocess.run()
    with (
        pytest.raises(SystemExit),
        patch("sys.argv", ["slurm-submit"]),
        patch("matbench_discovery.hpc.subprocess.run") as mock_subprocess_run,
        patch(
            "matbench_discovery.hpc.tempfile.mkdtemp",
            return_value="/tmp/slurm_job_123",
        ),
        patch("matbench_discovery.hpc.shutil.copy2") as mock_copy2,
        patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True),
    ):
        hpc.slurm_submit(**kwargs)  # ty: ignore[invalid-argument-type]
    assert mock_subprocess_run.call_count == 1

    expected_py_file_path = py_file_path or __file__
    if submit_as_temp_file:
        temp_path = f"/tmp/slurm_job_123/{os.path.basename(expected_py_file_path)}"
        mock_copy2.assert_called_once_with(expected_py_file_path, temp_path)
        expected_py_file_path = temp_path
    else:
        mock_copy2.assert_not_called()

    log_name = "slurm-%A-%a.log" if array else "slurm-%A.log"
    wrap_cmd = f"{expected_pre_cmd} python {expected_py_file_path}".strip()
    sbatch_cmd = (
        f"sbatch --job-name {job_name} --output {out_dir}/{log_name} --foo "
        f"--wrap {wrap_cmd}"
    ).replace(" --", "\n  --")
    for key, val in (
        ("time", time),
        ("account", account),
        ("partition", partition),
        ("array", array),
    ):
        if val is not None:
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
    ("segments", "expected"),
    [
        (
            [
                # unrelated keys like excluded_formula_reasons are ignored by merge
                {
                    "hardware": "H200",
                    "run_time_sec": 100,
                    "excluded_formula_reasons": {},
                },
                {"hardware": "H200", "run_time_sec": 50.5},
            ],
            {"hardware": "H200", "run_time_sec": 150.5},
        ),
        ([], {}),
        ([{}, {}], {}),  # no provenance -> empty (leaves existing YAML untouched)
        # all-or-nothing: a field missing from any segment is dropped, not summed
        # over the segments that have it
        ([{"run_time_sec": 12.0}, {}], {}),
        (
            [
                {"hardware": "H200", "run_time_sec": 10.0, "max_rss_gb": 4.0},
                {"hardware": "H200", "run_time_sec": 20.0},
            ],
            {"hardware": "H200", "run_time_sec": 30.0},
        ),
        (
            # memory peaks max over segments (each segment is its own process)
            [
                {"run_time_sec": 10.0, "max_rss_gb": 4.2, "max_gpu_mem_gb": 11.5},
                {"run_time_sec": 20.0, "max_rss_gb": 6.1, "max_gpu_mem_gb": 9.0},
            ],
            {"run_time_sec": 30.0, "max_rss_gb": 6.1, "max_gpu_mem_gb": 11.5},
        ),
        # mixed or missing GPU labels drop hardware instead of keeping the first one
        ([{"hardware": "H200"}, {"hardware": "H100"}], {}),
        ([{"hardware": "H200"}, {"hardware": ""}], {}),
        # NaN counts as missing, so the result can't depend on segment order
        ([{"max_rss_gb": float("nan")}, {"max_rss_gb": 2.0}], {}),
        ([{"max_rss_gb": 2.0}, {"max_rss_gb": float("nan")}], {}),
        # versions only when every segment agrees
        ([{"versions": {"numpy": "1"}}, {}], {}),
        ([{"versions": {"numpy": "1"}}, {"versions": {"numpy": "2"}}], {}),
        (
            [{"versions": {"numpy": "1"}}, {"versions": {"numpy": "1"}}],
            {"versions": {"numpy": "1"}},
        ),
        (
            # latest completion time and sorted, de-duplicated host/job unions
            [
                {
                    "completed_at": "2026-07-02T00:00:00",
                    "hostname": "node-b",
                    "slurm_job_id": "2",
                    "slurm_array_task_id": "1",
                },
                {
                    "completed_at": "2026-07-03T00:00:00",
                    "hostname": "node-a",
                    "slurm_job_id": "2",
                    "slurm_array_task_id": "0",
                },
                {"completed_at": "2026-07-01T00:00:00", "hostname": "node-b"},
            ],
            {
                "completed_at": "2026-07-03T00:00:00",
                "hostnames": ["node-a", "node-b"],
                "slurm_job_ids": ["2"],
                "slurm_array_task_ids": ["0", "1"],
            },
        ),
    ],
)
def test_merge_audit_metadata(
    segments: list[dict[str, object]], expected: dict[str, object]
) -> None:
    """Cost fields merge all-or-nothing; audit fields take latest time + unions."""
    assert hpc.merge_audit_metadata(segments) == expected


STRICT_SEGMENT = {"hardware": "H200", "run_time_sec": 1.0, "versions": {"numpy": "1"}}


@pytest.mark.parametrize(
    ("other_segment", "error_match"),
    [
        ({**STRICT_SEGMENT, "versions": {"numpy": "2"}}, "mixed or missing package"),
        ({**STRICT_SEGMENT, "hardware": ""}, "mixed or missing hardware"),
        ({**STRICT_SEGMENT, "hardware": "H100"}, "mixed or missing hardware"),
        ({**STRICT_SEGMENT, "run_time_sec": None}, "missing runtime"),
        ({**STRICT_SEGMENT, "max_rss_gb": 3.0}, "partial max_rss_gb"),
    ],
)
def test_merge_audit_metadata_strict_errors(
    other_segment: dict[str, object], error_match: str
) -> None:
    """Strict merging rejects mixed or incomplete provenance."""
    assert hpc.merge_audit_metadata([STRICT_SEGMENT] * 2, strict=True) == {
        "hardware": "H200",
        "run_time_sec": 2.0,
        "versions": {"numpy": "1"},
    }
    with pytest.raises(ValueError, match=error_match):
        hpc.merge_audit_metadata([STRICT_SEGMENT, other_segment], strict=True)


def test_peak_memory_gb() -> None:
    """peak_memory_gb reports a positive Unix RSS peak; GPU reset is a safe no-op."""
    mem = hpc.peak_memory_gb()
    if sys.platform != "win32":  # Windows lacks the resource module
        assert mem["max_rss_gb"] > 0
    assert set(mem) <= {"max_rss_gb", "max_gpu_mem_gb"}
    hpc.reset_gpu_peak_memory()  # must not raise without torch/CUDA
