import sys
from unittest.mock import patch

import pytest
from pytest import CaptureFixture

from mb_discovery.slurm import _get_calling_file_path, slurm_submit_python


def test_slurm_submit(capsys: CaptureFixture[str]) -> None:
    # check slurm_submit_python() does nothing in normal mode

    kwargs = dict(
        job_name="test_job",
        log_dir="test_log_dir",
        time="0:0:1",
        slurm_flags=("--test-flag",),
    )
    slurm_submit_python(**kwargs)  # type: ignore

    stdout, stderr = capsys.readouterr()
    assert stderr == stderr == ""

    # check slurm_submit_python() does prints cmd calls subprocess.run() in submit mode
    sys.argv += ["slurm-submit"]

    with pytest.raises(SystemExit), patch(
        "mb_discovery.slurm.subprocess.run"
    ) as mock_subprocess_run:
        slurm_submit_python(**kwargs)  # type: ignore

    assert mock_subprocess_run.call_count == 1

    sbatch_cmd = (
        "sbatch --partition=icelake --account=LEE-SL3-CPU --time=0:0:1 --job-name "
        "test_job --output test_log_dir/slurm-%A-%a.out --test-flag --wrap python"
    )
    stdout, stderr = capsys.readouterr()
    assert sbatch_cmd in stdout
    assert stderr == ""


def test_get_calling_file_path() -> None:
    assert _get_calling_file_path(frame=1) == __file__

    def wrapper(frame: int) -> str:
        return _get_calling_file_path(frame)

    assert wrapper(frame=2) == __file__
