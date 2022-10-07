import os
import subprocess
import sys
from collections.abc import Sequence


def _get_calling_file_path(frame: int = 1) -> str:
    """Return calling file's path.

    Args:
        frame (int, optional): How many function call's up? Defaults to 1.

    Returns:
        str: Calling function's file path n frames up the stack.
    """
    caller_path = sys._getframe(frame).f_code.co_filename
    return os.path.abspath(caller_path)


def slurm_submit_python(
    job_name: str,
    log_dir: str,
    time: str,
    py_file_path: str = None,
    slurm_flags: Sequence[str] = (),
    partition: str = "icelake",
    account: str = "LEE-SL3-CPU",
    array: str = "",
    pre_cmd: str = "",
) -> None:
    """Slurm submit a python script using sbatch --wrap 'python path/to/file.py' by
    calling this function in the script and invoking the script with
    `python path/to/file.py slurm-submit`.

    Args:
        job_name (str): Slurm job name.
        log_dir (str): Directory to write slurm logs. Log file will include job ID and
            array task ID.
        time (str): 'HH:MM:SS' time limit for the job.
        py_file_path (str): Path to the python script to be submitted. Defaults to the
            path of the file calling slurm_submit_python().
        slurm_flags (Sequence[str], optional): Extra slurm CLI flags. Defaults to ().
        partition (str, optional): Slurm partition. Defaults to "icelake".
        account (str, optional): Account to charge for this job.
            Defaults to "LEE-SL3-CPU".
        array (str, optional): Slurm array specifier. Defaults to "".
        pre_cmd (str, optional): Things like `module load` commands and environment
            variables to set when running the python script go here. Example:
            pre_cmd='ENV_VAR=42' or 'module load rhel8/default-amp;'. Defaults to "".

    Raises:
        SystemExit: Exit code will be subprocess.run(['sbatch', ...]).returncode.
    """
    if py_file_path is None:
        py_file_path = _get_calling_file_path(frame=2)

    cmd = [
        *f"sbatch --{partition=} --{account=} --{time=}".replace("'", "").split(),
        *("--job-name", job_name),
        *("--output", f"{log_dir}/slurm-%A-%a.out"),
        *slurm_flags,
        *("--wrap", f"{pre_cmd} python {py_file_path}".strip()),
    ]
    if array:
        cmd += ["--array", array]

    running_as_slurm_job = "SLURM_JOB_ID" in os.environ
    if running_as_slurm_job or "slurm-submit" in sys.argv:
        # print sbatch command at submission time and into slurm log file
        # but not when running interactively
        print(" ".join(cmd))

    if "slurm-submit" not in sys.argv:
        return

    os.makedirs(log_dir, exist_ok=True)  # slurm fails if log_dir is missing

    result = subprocess.run(cmd, check=True)

    raise SystemExit(result.returncode)
