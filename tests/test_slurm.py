from __future__ import annotations

import os
from unittest.mock import patch

import pytest
from pytest import CaptureFixture

from matbench_discovery.slurm import _get_calling_file_path, slurm_submit


@patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True)
@pytest.mark.parametrize("py_file_path", [None, "path/to/file.py"])
def test_slurm_submit(capsys: CaptureFixture[str], py_file_path: str | None) -> None:
    job_name = "test_job"
    out_dir = "tmp"
    time = "0:0:1"
    partition = "fake-partition"
    account = "fake-account"

    func_call = lambda: slurm_submit(
        job_name=job_name,
        out_dir=out_dir,
        time=time,
        partition=partition,
        account=account,
        py_file_path=py_file_path,
        slurm_flags=("--test-flag",),
    )

    slurm_vars = func_call()

    assert slurm_vars == {"slurm_job_id": "1234"}
    stdout, stderr = capsys.readouterr()
    # check slurm_submit() did nothing in normal mode
    assert stderr == stderr == ""

    # check slurm_submit() prints cmd and calls subprocess.run() in submit mode
    with pytest.raises(SystemExit), patch("sys.argv", ["slurm-submit"]), patch(
        "matbench_discovery.slurm.subprocess.run"
    ) as mock_subprocess_run:
        func_call()

    assert mock_subprocess_run.call_count == 1

    sbatch_cmd = (
        f"sbatch --partition={partition} --account={account} --time={time} "
        f"--job-name {job_name} --output {out_dir}/slurm-%A.log --test-flag "
        f"--wrap python {py_file_path or __file__}"
    ).replace(" --", "\n  --")
    stdout, stderr = capsys.readouterr()
    assert sbatch_cmd in stdout
    assert stderr == ""


def test_get_calling_file_path() -> None:
    assert _get_calling_file_path(frame=1) == __file__

    def wrapper(frame: int) -> str:
        return _get_calling_file_path(frame)

    assert wrapper(frame=2) == __file__
