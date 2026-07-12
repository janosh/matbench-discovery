"""Slurm job submission helper function."""

import importlib.metadata
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Mapping, Sequence, Sized
from typing import TYPE_CHECKING, Final, TypeVar

import numpy as np

if TYPE_CHECKING:
    import pandas as pd


def detect_hardware() -> str:
    """Human-readable name for the accelerator the run executed on.

    Tries torch, then JAX, then TensorFlow (MLIP backends expose the GPU differently),
    falling back to the CPU model when no GPU is visible or accelerator probes fail.
    """
    try:
        import torch

        if torch.cuda.is_available():
            return torch.cuda.get_device_name(0)
    except Exception:  # noqa: BLE001, S110
        pass
    try:
        import jax

        if gpus := [dev for dev in jax.devices() if dev.platform == "gpu"]:
            return gpus[0].device_kind
    except Exception:  # noqa: BLE001, S110
        pass
    try:
        import tensorflow as tf

        if gpus := tf.config.list_physical_devices("GPU"):
            details = tf.config.experimental.get_device_details(gpus[0])
            return details.get("device_name", "GPU")
    except Exception:  # noqa: BLE001, S110
        pass
    return f"CPU ({platform.processor() or platform.machine()})"


def reset_gpu_peak_memory() -> None:
    """Reset torch's peak-VRAM counter so the next peak_memory_gb() call attributes
    GPU memory to work done after this point (e.g. one system's rollout in a
    multi-system process). No-op without torch/CUDA; host RSS has no reset (it is a
    process-lifetime high-water mark).
    """
    # broad except: a broken CUDA runtime (e.g. driver mismatch mid-job) raises
    # RuntimeError, which must not crash the rollout this instruments
    try:
        import torch

        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
    except Exception:  # noqa: BLE001, S110
        pass


def peak_memory_gb() -> dict[str, float]:
    """Peak memory footprint of this process in GB (1e9 bytes), rounded to 3 decimals.

    Returns up to two entries: ``max_rss_gb``, the host resident-set-size high-water
    mark since process start (via getrusage; omitted on platforms without the
    ``resource`` module, i.e. Windows), and ``max_gpu_mem_gb``, torch's peak allocated
    CUDA memory since the last reset_gpu_peak_memory() call (omitted without
    torch/CUDA).
    """
    mem: dict[str, float] = {}
    try:
        import resource

        # ru_maxrss unit is bytes on macOS but KiB on Linux
        factor = 1 if sys.platform == "darwin" else 1024
        max_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * factor
        mem["max_rss_gb"] = round(max_rss / 1e9, 3)
    except ImportError:
        pass
    # broad except like reset_gpu_peak_memory: memory introspection failing (broken
    # CUDA runtime) must not crash collect_run_info after a multi-hour rollout
    try:
        import torch

        if torch.cuda.is_available():
            mem["max_gpu_mem_gb"] = round(torch.cuda.max_memory_allocated() / 1e9, 3)
    except Exception:  # noqa: BLE001, S110
        pass
    return mem


def package_versions(package_names: Iterable[str]) -> dict[str, str]:
    """Return installed versions for Python and the requested packages."""
    versions = {"python": sys.version.split()[0]}
    for package_name in sorted(set(package_names)):
        try:
            versions[package_name] = importlib.metadata.version(package_name)
        except importlib.metadata.PackageNotFoundError:
            continue
    return versions


# taken from https://slurm.schedmd.com/job_array.html#env_vars, lower-cased and
# and removed the SLURM_ prefix
SLURM_KEYS: Final[tuple[str, ...]] = (
    "job_id",
    "array_job_id",
    "array_task_id",
    "array_task_count",
    "mem_per_node",
    "nodelistsubmit_host",
    "job_partition",
    "job_user",
    "job_account",
    "tasks_per_node",
    "job_qos",
)
SLURM_SUBMIT_KEY: Final[str] = "slurm-submit"
HasLen = TypeVar("HasLen", bound=Sized)
COST_PROVENANCE_KEYS = ("hardware", "run_time_sec", "max_rss_gb", "max_gpu_mem_gb")


def slurm_shard_selection(
    n_shards_arg: int | None, shard_index_arg: int | None
) -> tuple[int, int]:
    """Resolve zero-based shard selection from CLI flags or Slurm environment."""
    n_shards = (
        n_shards_arg
        if n_shards_arg is not None
        else int(os.getenv("SLURM_ARRAY_TASK_COUNT") or 1)
    )
    if n_shards < 1:
        raise ValueError(f"n_shards must be positive, got {n_shards}")

    if shard_index_arg is not None:
        shard_index = shard_index_arg
    elif slurm_task_id := os.getenv("SLURM_ARRAY_TASK_ID"):
        if n_shards_arg is not None:
            # Explicit shard counts make partial reruns such as --array=3,7 map back
            # to their original zero-based shard indices.
            slurm_task_max = os.getenv("SLURM_ARRAY_TASK_MAX")
            if slurm_task_max is not None and int(slurm_task_max) >= n_shards:
                # fail every task of a 1-based --array=1-N immediately instead of
                # running N-1 shards and silently never producing shard 0
                raise ValueError(
                    f"SLURM_ARRAY_TASK_MAX={slurm_task_max} exceeds the last shard "
                    f"index {n_shards - 1}; with explicit --n-shards, array task IDs "
                    "are treated as zero-based original shard indices"
                )
            shard_index = int(slurm_task_id)
        else:
            slurm_task_min = int(os.getenv("SLURM_ARRAY_TASK_MIN", "0"))
            slurm_task_max_id = int(
                os.getenv("SLURM_ARRAY_TASK_MAX", str(slurm_task_min + n_shards - 1))
            )
            if slurm_task_max_id - slurm_task_min + 1 != n_shards:
                raise ValueError(
                    "Partial Slurm arrays require explicit --n-shards so task IDs "
                    "retain their original shard indices (any rerun of a shard "
                    "subset, contiguous or not, must pass the original --n-shards)"
                )
            shard_index = int(slurm_task_id) - slurm_task_min
    elif n_shards == 1:
        shard_index = 0
    else:
        raise ValueError("--shard-index is required outside a Slurm array")

    if not 0 <= shard_index < n_shards:
        raise ValueError(
            f"shard_index must be in [0, {n_shards - 1}], got {shard_index}"
        )
    return n_shards, shard_index


def effective_shard_args(
    n_shards: int | None,
    shard_index: int | None,
    *,
    dry_run: bool,
) -> tuple[int | None, int | None, bool]:
    """Run one dry-run shard on only the first task of an implicit Slurm array."""
    if (
        not dry_run
        or n_shards is not None
        or shard_index is not None
        or not os.getenv("SLURM_ARRAY_TASK_COUNT")
    ):
        return n_shards, shard_index, False
    slurm_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
    slurm_task_min = int(os.getenv("SLURM_ARRAY_TASK_MIN", "0"))
    return 1, 0, slurm_task_id != slurm_task_min


def partition_material_ids(
    items_by_id: Mapping[str, Sized], n_shards: int
) -> list[list[str]]:
    """Deterministically balance sized items across shards by total length."""
    if n_shards < 1:
        raise ValueError(f"n_shards must be positive, got {n_shards}")
    if n_shards > len(items_by_id):
        raise ValueError(
            f"n_shards ({n_shards}) cannot exceed structures ({len(items_by_id)})"
        )
    shards: list[list[str]] = [[] for _ in range(n_shards)]
    shard_sizes = [0] * n_shards
    for item_id, item in sorted(
        items_by_id.items(), key=lambda pair: (-len(pair[1]), pair[0])
    ):
        shard_idx = min(range(n_shards), key=lambda idx: (shard_sizes[idx], idx))
        shards[shard_idx].append(item_id)
        shard_sizes[shard_idx] += len(item)
    return shards


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


def merge_run_metadata(
    shard_metadatas: list[dict[str, object]],
) -> dict[str, str | float]:
    """Shared hardware, summed run_time_sec and peak memory across Slurm shard
    run_metadata dicts.

    Shards run in parallel, so the summed run_time_sec ~ a serial sweep, memory peaks
    are maxed (each shard is its own process) and hardware is shared. Keys no shard
    recorded are omitted, so a merge keeps existing YAML values.
    """
    hardware = next(
        (meta["hardware"] for meta in shard_metadatas if meta.get("hardware")), None
    )
    run_times = [
        float(secs)
        for meta in shard_metadatas
        if isinstance(secs := meta.get("run_time_sec"), int | float)
    ]
    merged: dict[str, str | float] = {
        **({"hardware": hardware} if isinstance(hardware, str) else {}),
        **({"run_time_sec": round(sum(run_times), 2)} if run_times else {}),
    }
    for mem_key in ("max_rss_gb", "max_gpu_mem_gb"):
        peaks = [
            val
            for meta in shard_metadatas
            if isinstance(val := meta.get(mem_key), int | float)
        ]
        if peaks:
            merged[mem_key] = max(peaks)
    return merged


def merge_audit_metadata(
    metadata_segments: Sequence[Mapping[str, object]], *, strict: bool = False
) -> dict[str, object]:
    """Merge cost and audit provenance, optionally rejecting incomplete segments."""
    segments = [dict(segment) for segment in metadata_segments]
    versions = [
        value
        for segment in segments
        if isinstance(value := segment.get("versions"), dict)
    ]
    hardware_labels = {
        str(segment["hardware"]) for segment in segments if segment.get("hardware")
    }
    cost_counts = {
        key: sum(
            (
                isinstance(segment.get(key), str) and bool(segment.get(key))
                if key == "hardware"
                else isinstance(segment.get(key), int | float)
            )
            for segment in segments
        )
        for key in COST_PROVENANCE_KEYS
    }
    versions_complete_and_equal = len(versions) == len(segments) and all(
        value == versions[0] for value in versions[1:]
    )
    if strict:
        if not versions_complete_and_equal:
            raise ValueError("Run segments contain mixed or missing package versions")
        if cost_counts["hardware"] != len(segments) or len(hardware_labels) != 1:
            raise ValueError(
                "Run segments contain mixed or missing hardware provenance"
            )
        if cost_counts["run_time_sec"] != len(segments):
            raise ValueError("Run segments contain missing runtime provenance")
        for memory_key in COST_PROVENANCE_KEYS[2:]:
            if cost_counts[memory_key] not in {0, len(segments)}:
                raise ValueError(
                    f"Run segments contain partial {memory_key} provenance"
                )

    merged: dict[str, object] = dict(merge_run_metadata(segments))
    if len(hardware_labels) > 1:
        merged.pop("hardware", None)
    for cost_key, count in cost_counts.items():
        if count != len(segments):
            merged.pop(cost_key, None)
    if versions and versions_complete_and_equal:
        merged["versions"] = versions[0]
    if completed_at := [
        str(segment["completed_at"])
        for segment in segments
        if segment.get("completed_at")
    ]:
        merged["completed_at"] = max(completed_at)
    for source_key, merged_key in (
        ("hostname", "hostnames"),
        ("slurm_job_id", "slurm_job_ids"),
        ("slurm_array_task_id", "slurm_array_task_ids"),
    ):
        values = sorted(
            {
                str(segment[source_key])
                for segment in segments
                if segment.get(source_key)
            }
        )
        if values:
            merged[merged_key] = values
    return merged


def df_slurm_chunk(
    df_in: "pd.DataFrame", n_chunks: int, task_id: int
) -> "pd.DataFrame":
    """Get the chunk a 1-based slurm array task should process, i.e. the
    (task_id - 1)-th of n_chunks roughly equal row-wise splits of df_in.

    np.array_split(df_in, ...) no longer works since pandas 3 removed
    DataFrame.swapaxes which numpy relied on to preserve the DataFrame type, so
    split row indices instead and select with iloc.
    """
    if n_chunks < 1 or not 1 <= task_id <= n_chunks:
        raise ValueError(f"expected {n_chunks=} >= 1 and 1 <= {task_id=} <= n_chunks")
    row_indices = np.array_split(np.arange(len(df_in)), n_chunks)[task_id - 1]
    return df_in.iloc[row_indices]


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
    chunk_sizes = np.zeros(n_chunks, dtype=float)

    # Assign each structure to the chunk with the smallest current total
    for sized_obj in sorted_inputs:
        smallest_chunk_idx = int(np.argmin(chunk_sizes))
        chunks[smallest_chunk_idx].append(sized_obj)
        chunk_sizes[smallest_chunk_idx] += len(sized_obj)

    if report:
        # Print statistics about the chunk sizes
        mean, std = chunk_sizes.mean(), chunk_sizes.std()
        min_chunk_size = float(np.min(chunk_sizes))
        max_chunk_size = float(np.max(chunk_sizes))
        cls_name = type(inputs[0]).__name__
        print(
            f"Split {len(inputs):,} structures into {n_chunks:,} chunks:\n"
            f"Mean sum(len({cls_name})) per chunk: {mean:,.1f} ± {std:,.1f}, "
            f"min: {min_chunk_size:,.0f}, max: {max_chunk_size:,.0f}"
        )

    return chunks
