from __future__ import annotations

from unittest.mock import patch

import pytest
from pytest import CaptureFixture

from mb_discovery.slurm import _get_calling_file_path, slurm_submit_python


@pytest.mark.parametrize("py_file_path", [None, "path/to/file.py"])
def test_slurm_submit(capsys: CaptureFixture[str], py_file_path: str | None) -> None:
    kwargs = dict(
        job_name="test_job",
        log_dir="test_log_dir",
        partition="fake-partition",
        account="fake-account",
        py_file_path=py_file_path,
        time="0:0:1",
        slurm_flags=("--test-flag",),
    )
    slurm_submit_python(**kwargs)  # type: ignore

    stdout, stderr = capsys.readouterr()
    # check slurm_submit_python() did nothing in normal mode
    assert stderr == stderr == ""

    # check slurm_submit_python() prints cmd and calls subprocess.run() in submit mode
    with pytest.raises(SystemExit), patch("sys.argv", ["slurm-submit"]), patch(
        "mb_discovery.slurm.subprocess.run"
    ) as mock_subprocess_run:
        slurm_submit_python(**kwargs)  # type: ignore

    assert mock_subprocess_run.call_count == 1

    sbatch_cmd = (
        "sbatch --partition=fake-partition --account=fake-account --time=0:0:1 "
        "--job-name test_job --output test_log_dir/slurm-%A-%a.out --test-flag "
        f"--wrap python {py_file_path or __file__}"
    )
    stdout, stderr = capsys.readouterr()
    assert sbatch_cmd in stdout
    assert stderr == ""


def test_get_calling_file_path() -> None:
    assert _get_calling_file_path(frame=1) == __file__

    def wrapper(frame: int) -> str:
        return _get_calling_file_path(frame)

    assert wrapper(frame=2) == __file__
