from __future__ import annotations

import os
import subprocess
import sys
from collections.abc import Sequence

SLURM_KEYS = (
    "job_id array_task_id array_task_count mem_per_node nodelist submit_host"
    "job_partition job_user job_account tasks_per_node job_qos"
).split()


def _get_calling_file_path(frame: int = 1) -> str:
    """Return calling file's path.

    Args:
        frame (int, optional): How many function call's up? Defaults to 1.

    Returns:
        str: Calling function's file path n frames up the stack.
    """
    caller_path = sys._getframe(frame).f_code.co_filename
    return os.path.abspath(caller_path)


def slurm_submit(
    job_name: str,
    out_dir: str,
    time: str,
    partition: str,
    account: str,
    py_file_path: str | None = None,
    slurm_flags: str | Sequence[str] = (),
    array: str | None = None,
    pre_cmd: str = "",
) -> dict[str, str]:
    """Slurm submits a python script using `sbatch --wrap 'python path/to/file.py'`.

    Usage: Call this function at the top of the script (before doing any real work) and
    then submit a job with `python path/to/that/script.py slurm-submit`. The slurm job
    will run the whole script.

    Args:
        job_name (str): Slurm job name.
        out_dir (str): Directory to write slurm logs. Log file will include slurm job
            ID and array task ID.
        time (str): 'HH:MM:SS' time limit for the job.
        py_file_path (str, optional): Path to the python script to be submitted.
            Defaults to the path of the file calling slurm_submit().
        partition (str, optional): Slurm partition.
        account (str, optional): Account to charge for this job.
        slurm_flags (str | list[str], optional): Extra slurm CLI flags. Defaults to ().
            Examples: ('--nodes 1', '--gpus-per-node 1') or ('--mem', '16G').
        array (str, optional): Slurm array specifier. Defaults to None. Example:
            '9' (for SLURM_ARRAY_TASK_ID from 0-9 inclusive), '1-10' or '1-10%2', etc.
        pre_cmd (str, optional): Things like `module load` commands and environment
            variables to set when running the python script go here. Example:
            pre_cmd='ENV_VAR=42' or 'module load rhel8/default-amp;'. Defaults to "".
            If running on CPU, pre_cmd="unset OMP_NUM_THREADS" allows PyTorch to use
            all cores https://docs.hpc.cam.ac.uk/hpc/software-packages/pytorch.html

    Raises:
        SystemExit: Exit code will be subprocess.run(['sbatch', ...]).returncode.

    Returns:
        dict[str, str]: Slurm variables like job ID, array task ID, compute nodes IDs,
            submission node ID and total job memory.
    """
    if py_file_path is None:
        py_file_path = _get_calling_file_path(frame=2)

    if "GPU" in partition:
        # on Ampere GPU partition, source module CLI and load default Ampere env
        # before actual job command
        pre_cmd += ". /etc/profile.d/modules.sh; module load rhel8/default-amp;"

    os.makedirs(out_dir, exist_ok=True)  # slurm fails if out_dir is missing

    cmd = [
        *f"sbatch --{partition=} --{account=} --{time=}".replace("'", "").split(),
        *("--job-name", job_name),
        *("--output", f"{out_dir}/slurm-%A{'-%a' if array else ''}.log"),
        *(slurm_flags.split() if isinstance(slurm_flags, str) else slurm_flags),
        *("--wrap", f"{pre_cmd} python {py_file_path}".strip()),
    ]
    if array:
        cmd += ["--array", array]

    is_log_file = not sys.stdout.isatty()
    is_slurm_job = "SLURM_JOB_ID" in os.environ

    slurm_vars = {
        f"slurm_{key}": val
        for key in SLURM_KEYS
        if (val := os.getenv(f"SLURM_{key}".upper()))
    }
    slurm_vars["slurm_timelimit"] = time
    if slurm_flags:
        slurm_vars["slurm_flags"] = str(slurm_flags)
    if pre_cmd:
        slurm_vars["pre_cmd"] = pre_cmd

    # print sbatch command into slurm log file and at job submission time
    # but not into terminal or Jupyter
    if (is_slurm_job and is_log_file) or "slurm-submit" in sys.argv:
        print(f"\n{' '.join(cmd)}\n".replace(" --", "\n  --"))
    if is_slurm_job and is_log_file:
        for key, val in slurm_vars.items():
            print(f"{key}={val}")

    if "slurm-submit" not in sys.argv:
        return slurm_vars  # if not submitting slurm job, resume outside code as normal

    result = subprocess.run(cmd, check=True)

    # after sbatch submission, exit with slurm exit code
    raise SystemExit(result.returncode)
