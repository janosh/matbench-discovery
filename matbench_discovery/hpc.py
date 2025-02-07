"""Slurm job submission helper function."""

import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Sequence, Sized
from typing import TypeVar

import numpy as np

# taken from https://slurm.schedmd.com/job_array.html#env_vars, lower-cased and
# and removed the SLURM_ prefix
SLURM_KEYS = (
    "job_id array_job_id array_task_id array_task_count mem_per_node nodelist"
    "submit_host job_partition job_user job_account tasks_per_node job_qos"
).split()
SLURM_SUBMIT_KEY = "slurm-submit"

HasLen = TypeVar("HasLen", bound=Sized)


def _get_calling_file_path(frame: int = 1) -> str:
    """Return calling file's path.

    Args:
        frame (int, optional): How many function call's up? Defaults to 1.

    Returns:
        str: Calling function's file path n frames up the stack.
    """
    caller_path = sys._getframe(frame).f_code.co_filename  # noqa: SLF001
    return os.path.abspath(caller_path)


def slurm_submit(
    job_name: str,
    out_dir: str,
    *,
    time: str | None = None,
    account: str | None = None,
    partition: str | None = None,
    py_file_path: str | None = None,
    slurm_flags: str | Sequence[str] = (),
    array: str | None = None,
    pre_cmd: str = "",
    submit_as_temp_file: bool = True,
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
            Defaults to the path of the file calling slurm_submit().
        account (str): Account to charge for this job.
        partition (str, optional): Slurm partition.
        py_file_path (str, optional): Path to the python script to be submitted.
        slurm_flags (str | list[str], optional): Extra slurm CLI flags. Defaults to ().
            Examples: ('--nodes 1', '--gpus-per-node 1') or ('--mem', '16G').
        array (str, optional): Slurm array specifier. Defaults to None. Example:
            '9' (for SLURM_ARRAY_TASK_ID from 0-9 inclusive), '1-10' or '1-10%2', etc.
        pre_cmd (str, optional): Things like `module load` commands and environment
            variables to set before running the python script go here. Example:
            pre_cmd='ENV_VAR=42' or 'module load pytorch;'. Defaults to "". If running
            on CPU, pre_cmd="unset OMP_NUM_THREADS" allows PyTorch to use all cores.
        submit_as_temp_file (bool, optional): If True, copy the Python file to a
            temporary directory before submitting. This allows the user to modify
            the original file without affecting queued jobs. Defaults to True.

    Raises:
        SystemExit: Exit code will be subprocess.run(['sbatch', ...]).returncode.

    Returns:
        dict[str, str]: Slurm variables like job ID, array task ID, compute nodes IDs,
            submission node ID and total job memory.
    """
    py_file_path = py_file_path or _get_calling_file_path(frame=2)

    os.makedirs(out_dir, exist_ok=True)  # slurm fails if out_dir is missing

    # Copy the file to a temporary directory if submit_as_temp_file is True
    if submit_as_temp_file and SLURM_SUBMIT_KEY in sys.argv:
        temp_dir = tempfile.mkdtemp(prefix="slurm_job_")
        temp_file_path = f"{temp_dir}/{os.path.basename(py_file_path)}"
        shutil.copy2(py_file_path, temp_file_path)
        py_file_path = temp_file_path

    # ensure pre_cmd ends with a semicolon
    if pre_cmd and not pre_cmd.strip().endswith(";"):
        pre_cmd += ";"

    cmd = [
        *("sbatch", "--job-name", job_name),
        *("--output", f"{out_dir}/slurm-%A{'-%a' if array else ''}.log"),
        *(slurm_flags.split() if isinstance(slurm_flags, str) else slurm_flags),
        *("--wrap", f"{pre_cmd or ''} python {py_file_path}".strip()),
    ]
    for flag in (f"{time=!s}", f"{account=!s}", f"{partition=!s}", f"{array=!s}"):
        key, val = flag.split("=")
        if val != "None":
            cmd += (f"--{key}", val)

    is_log_file = not sys.stdout.isatty()
    is_slurm_job = "SLURM_JOB_ID" in os.environ

    slurm_vars = {
        f"slurm_{key}": os.environ[f"SLURM_{key}".upper()]
        for key in SLURM_KEYS
        if f"SLURM_{key}".upper() in os.environ
    }
    if time is not None:
        slurm_vars["slurm_timelimit"] = time
    if slurm_flags != ():
        slurm_vars["slurm_flags"] = str(slurm_flags)
    if pre_cmd not in ("", None):
        slurm_vars["pre_cmd"] = pre_cmd

    # print sbatch command into slurm log file and at job submission time
    # but not into terminal or Jupyter
    if (is_slurm_job and is_log_file) or SLURM_SUBMIT_KEY in sys.argv:
        print(f"\n{' '.join(cmd)}\n".replace(" --", "\n  --"))
    if is_slurm_job and is_log_file:
        for key, val in slurm_vars.items():
            print(f"{key}={val}")

    if SLURM_SUBMIT_KEY not in sys.argv:
        return slurm_vars  # if not submitting slurm job, resume outside code as normal

    result = subprocess.run(cmd, check=True)

    # after sbatch submission, exit with slurm exit code
    raise SystemExit(result.returncode)


def chunk_by_lens(
    inputs: Sequence[HasLen],
    *,  # force keyword-only arguments
    n_chunks: int | None = None,
    chunk_size: int | None = None,
    report: bool = True,
) -> list[list[HasLen]]:
    """Make a balanced partition. That is, split a list of pymatgen Structures or
    ASE Atoms or anything with len() into chunks with roughly equal total length.

    This is useful for distributing workload evenly among workers, since computational
    cost often scales with the number of atoms in a structure.

    Args:
        inputs (Sequence[T]): List of objects with len() to split into chunks
        n_chunks (int, optional): Number of chunks to create. Defaults to None.
        chunk_size (int, optional): Target size for each chunk. Defaults to None.
            Only one of n_chunks or chunk_size can be specified.
        report (bool, optional): If True, print statistics about the chunk sizes.

    Returns:
        list[list[T]]: Each sublist contains objects and the sum of len() of objects
            in each chunk is roughly equal.

    Example:
        >>> from ase.build import bulk
        >>> structures = [bulk("Cu") * (i, i, 1) for i in range(1, 5)]
        >>> # Split into 2 chunks
        >>> chunks = chunk_by_lens(structures, n_chunks=2)
        >>> [sum(len(atoms) for atoms in chunk) for chunk in chunks]
        [14, 16]  # roughly equal total atom counts
        >>> # Or split into chunks of ~10 atoms each
        >>> chunks = chunk_by_lens(structures, chunk_size=10)
        >>> [sum(len(atoms) for atoms in chunk) for chunk in chunks]
        [12, 10, 8]  # roughly equal total atom counts

    Raises:
        ValueError: If neither or both n_chunks and chunk_size are specified.
    """
    if len(inputs) == 0:
        return []

    if n_chunks is not None and chunk_size is not None:
        raise ValueError("Cannot specify both n_chunks and chunk_size")

    # Get number of atoms in each structure
    lens = np.array([len(obj) for obj in inputs])
    total_size = lens.sum()

    if chunk_size:
        # Calculate n_chunks based on chunk_size
        n_chunks = max(1, int(np.ceil(total_size / chunk_size)))
    elif n_chunks:
        # n_chunks is specified
        if n_chunks < 1:
            raise ValueError("n_chunks must be >= 1")
        n_chunks = min(n_chunks, len(inputs))
    else:
        raise ValueError("n_chunks or chunk_size must be positive integer")

    # Sort structures by size (largest first) to help achieve better balance
    sort_idx = np.argsort(lens)[::-1]
    sorted_inputs = [inputs[i] for i in sort_idx]

    chunks: list[list[HasLen]] = [[] for _ in range(n_chunks)]
    chunk_sizes = np.zeros(n_chunks)

    # Assign each structure to the chunk with the smallest current total
    for sized_obj in sorted_inputs:
        smallest_chunk = np.argmin(chunk_sizes)
        chunks[smallest_chunk].append(sized_obj)
        chunk_sizes[smallest_chunk] += len(sized_obj)

    if report:
        # Print statistics about the chunk sizes
        mean, std = chunk_sizes.mean(), chunk_sizes.std()
        cls_name = type(inputs[0]).__name__
        print(
            f"Split {len(inputs):,} structures into {n_chunks:,} chunks:\n"
            f"Mean sum(len({cls_name})) per chunk: {mean:,.1f} ± {std:,.1f}, "
            f"min: {chunk_sizes.min():,.0f}, max: {chunk_sizes.max():,.0f}"
        )

    return chunks
