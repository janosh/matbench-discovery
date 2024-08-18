import os
from unittest.mock import patch

import pytest

from matbench_discovery.slurm import _get_calling_file_path, slurm_submit


@pytest.mark.parametrize("py_file_path", [None, "path/to/file.py"])
@pytest.mark.parametrize("partition", [None, "fake-partition"])
@pytest.mark.parametrize("time", [None, "0:0:1"])
@pytest.mark.parametrize("account", [None, "fake-account"])
@pytest.mark.parametrize("pre_cmd", [None, "module load pytorch;", "ENV_VAR=42"])
def test_slurm_submit(
    capsys: pytest.CaptureFixture[str],
    py_file_path: str | None,
    partition: str | None,
    time: str | None,
    account: str | None,
    pre_cmd: str | None,
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
    )

    slurm_submit(**kwargs)  # type: ignore[arg-type]

    stdout, stderr = capsys.readouterr()
    # check slurm_submit() did nothing in normal mode
    assert stdout == stderr == ""

    with patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True):
        slurm_vars = slurm_submit(**kwargs)  # type: ignore[arg-type]
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
        patch("matbench_discovery.slurm.subprocess.run") as mock_subprocess_run,
        patch.dict(os.environ, {"SLURM_JOB_ID": "1234"}, clear=True),
    ):
        slurm_submit(**kwargs)  # type: ignore[arg-type]

    assert mock_subprocess_run.call_count == 1
    sbatch_cmd = (
        f"sbatch --job-name {job_name} --output {out_dir}/slurm-%A.log --foo "
        f"--wrap {pre_cmd + ' ' if pre_cmd else ''}python {py_file_path or __file__}"
    ).replace(" --", "\n  --")
    for flag in (f"{time=!s}", f"{account=!s}", f"{partition=!s}"):
        key, val = flag.split("=")
        if val != "None":
            sbatch_cmd += f"\n  --{key} {val}"

    stdout, stderr = capsys.readouterr()
    assert sbatch_cmd in stdout
    assert stderr == ""


def test_get_calling_file_path() -> None:
    assert _get_calling_file_path(frame=1) == __file__

    def wrapper(frame: int) -> str:
        return _get_calling_file_path(frame)

    assert wrapper(frame=2) == __file__
