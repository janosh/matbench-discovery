"""HPC helpers: Slurm submission, sharding and run-cost provenance."""

import hashlib
import importlib.metadata
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import types
import warnings
from collections.abc import Iterable, Mapping, Sequence, Sized
from typing import TYPE_CHECKING, Any, Final, TypeIs

from pymatviz.enums import Key

if TYPE_CHECKING:
    from ase import Atoms

# === Run-cost provenance ===
# Where a cost block came from: measured by our current runners and recorded in the
# artifact (runner), our earlier runs predating that provenance, e.g. Slurm-log
# salvage or load-inclusive times (salvaged), or the model authors
COST_SOURCES: Final = ("runner", "salvaged", "author_reported")
# tags valid without artifact-recorded cost (YAML-only numbers)
YAML_ONLY_COST_SOURCES: Final = ("salvaged", "author_reported")
# where the timed work ran: merged units must agree on each, and resumed work must
# match the sessions before it (plus their package versions)
RUN_LABEL_KEYS: Final = ("hardware", "device", "dtype", "acceleration", "backend")
# the YAML `cost: {...}` mapping next to a task artifact, in canonical order. Flat
# artifact run metadata spells its source key cost_source
COST_KEYS: Final = (
    *("source", "hardware", "device", "dtype", "acceleration", "n_shards"),
    *("n_processes", "run_time_sec", "load_time_sec", "max_rss_gb", "max_gpu_mem_gb"),
)
# trustworthy keys of artifacts predating cost_source: older runners recorded the
# requested (not the actual) dtype/device
LEGACY_ARTIFACT_COST_KEYS: Final = (
    *("hardware", "n_shards", "run_time_sec", "max_rss_gb", "max_gpu_mem_gb"),
)
# only published when hardware is uniform + complete
GPU_DEPENDENT_COST_KEYS: Final = ("run_time_sec", "load_time_sec", "max_gpu_mem_gb")
# runtime overrides once a process's model load is charged to one of its units
# (rollout, element file, material): merges sum load time and process counts
LOAD_CHARGED: Final = {"load_time_sec": 0.0, "n_processes": 0}


def detect_hardware(device: str) -> str:
    """Name of the torch-style ``device`` (e.g. ``cuda:1``) as
    torch.cuda.get_device_name spells it, ``CPU (<arch>)`` for CPUs.

    Raises:
        ValueError: If a CUDA device is requested but torch sees no GPU.
    """
    if not device.startswith("cuda"):
        return f"CPU ({platform.processor() or platform.machine()})"
    import torch

    if not torch.cuda.is_available():
        raise ValueError(f"{device=} requested but torch.cuda.is_available() is False")
    return torch.cuda.get_device_name(torch.device(device))


def _torch_runtime(obj: object) -> dict[str, str] | None:
    """Device, dtype and hardware of the first floating-point torch parameter or
    buffer reachable from ``obj`` (breadth-first over attributes and container
    items), None when no torch module is reachable.
    """
    if "torch" not in sys.modules:  # never import torch for JAX/TensorFlow models
        return None
    import torch

    skip_types = (str, bytes, int, float, type, types.ModuleType, types.FunctionType)
    seen: set[int] = set()
    frontier: list[object] = [obj]
    for _depth in range(5):  # obj plus up to 4 attribute/container hops
        next_frontier: list[object] = []
        for item in frontier:
            if id(item) in seen or isinstance(item, skip_types) or len(seen) > 20_000:
                continue
            seen.add(id(item))
            if isinstance(item, torch.nn.Module):
                for tensor in (*item.parameters(), *item.buffers()):
                    if tensor.is_floating_point():
                        return {
                            "backend": "torch",
                            "device": tensor.device.type,
                            "torch_device": str(tensor.device),
                            "dtype": str(tensor.dtype).removeprefix("torch."),
                            "hardware": detect_hardware(str(tensor.device)),
                        }
            elif isinstance(item, dict):
                next_frontier.extend(item.values())
            elif isinstance(item, list | tuple):
                next_frontier.extend(item)
            elif hasattr(item, "__dict__") and not isinstance(item, torch.Tensor):
                next_frontier.extend(vars(item).values())
        frontier = next_frontier
    return None


def detect_runtime(calculator: object, *, requested_device: str) -> dict[str, str]:
    """Backend, device, hardware and (when known) dtype and ``torch_device`` a
    calculator actually runs on.

    Torch models are located through their parameters (so a factory ignoring the
    requested device is labeled truthfully), JAX/TensorFlow models through their
    framework's first visible GPU. Unreachable torch modules fall back to backend
    ``torch_probe`` labeled by the requested device (memory_peaks later requires
    measured VRAM to back a GPU label).
    """
    if (torch_runtime := _torch_runtime(calculator)) is not None:
        return torch_runtime
    backend = next((name for name in ("jax", "tensorflow") if name in sys.modules), "")
    if backend == "jax":
        import jax

        if gpus := [dev for dev in jax.devices() if dev.platform == "gpu"]:
            return {"backend": "jax", "device": "gpu", "hardware": gpus[0].device_kind}
    elif backend == "tensorflow":
        import tensorflow as tf

        if gpus := tf.config.list_physical_devices("GPU"):
            details = tf.config.experimental.get_device_details(gpus[0])
            if not (name := details.get("device_name")):
                raise ValueError(f"TensorFlow GPU {gpus[0]} reports no device_name")
            return {"backend": backend, "device": "gpu", "hardware": name}
    elif requested_device.startswith("cuda"):
        return {
            "backend": "torch_probe",
            "device": "cuda",
            "torch_device": requested_device,
            "hardware": detect_hardware(requested_device),
        }
    return {
        "backend": backend or "none",
        "device": "cpu",
        "hardware": detect_hardware("cpu"),
    }


def _is_torch_cuda_run(runtime: Mapping[str, Any]) -> bool:
    """Whether torch's CUDA caching-allocator counters measure this run's VRAM."""
    return runtime.get("backend") in ("torch", "torch_probe") and str(
        runtime.get("device", "")
    ).startswith("cuda")


def reset_gpu_peak_memory(runtime: Mapping[str, Any]) -> None:
    """Reset a torch CUDA run's peak-VRAM counter so memory_peaks() measures only
    work after this point (e.g. one rollout of a multi-system process).
    """
    if not _is_torch_cuda_run(runtime):
        return
    # broad except: a broken CUDA runtime (e.g. driver mismatch mid-job) raises
    # RuntimeError, which must not crash the rollout this instruments
    try:
        import torch

        torch.cuda.reset_peak_memory_stats(runtime.get("torch_device"))
    except Exception:  # noqa: BLE001, S110
        pass


def memory_peaks(runtime: Mapping[str, Any], *, succeeded: bool) -> dict[str, float]:
    """Peak memory in GB (1e9 bytes, 3 decimals): ``max_rss_gb``, the host RSS
    high-water mark since process start (not on Windows), and for torch CUDA runs
    ``max_gpu_mem_gb``, the peak allocated VRAM on ``torch_device`` since the last
    reset_gpu_peak_memory() (torch's counter reads 0 for JAX/TensorFlow).

    Raises:
        ValueError: If ``succeeded`` (a failed unit may never touch the GPU) but a
            torch CUDA run allocated no VRAM: the model is not on its labeled GPU.
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
    if not _is_torch_cuda_run(runtime):
        return mem
    # broad except like reset_gpu_peak_memory: memory introspection failing (broken
    # CUDA runtime) must not crash a run after a multi-hour rollout
    try:
        import torch

        if torch.cuda.is_available():
            peak = torch.cuda.max_memory_allocated(runtime.get("torch_device"))
            mem["max_gpu_mem_gb"] = round(peak / 1e9, 3)
    except Exception:  # noqa: BLE001, S110
        pass
    if succeeded and not mem.get("max_gpu_mem_gb", 0) > 0:
        raise ValueError(
            f"Run labeled {runtime.get('hardware')!r} on device "
            f"{runtime.get('device')!r} ({runtime.get('backend')} backend) allocated "
            f"no CUDA memory ({mem}) despite successful work: the model is not "
            "running on that GPU"
        )
    return mem


def check_resume_labels(
    previous: Mapping[str, Any], runtime: Mapping[str, Any], context: str
) -> None:
    """Require work resumed by ``runtime`` to share the run labels and package
    versions of the session that wrote ``previous``: a summed cost over two GPUs,
    precisions or environments describes neither.

    Raises:
        ValueError: Naming the differing labels (before, now).
    """
    keys = (*RUN_LABEL_KEYS, "versions")
    labels = [{key: meta.get(key) for key in keys} for meta in (previous, runtime)]
    # JSON-normalized: stored labels round-tripped through JSON (tuples -> lists)
    before, now = json.loads(json.dumps(labels, default=str))
    if diff := {
        key: (before[key], now[key]) for key in keys if before[key] != now[key]
    }:
        raise ValueError(
            f"{context}: cannot resume under different run labels (before, now): "
            f"{diff}; resume on the original hardware/environment or delete the "
            "partial outputs to restart"
        )


def package_versions(package_names: Iterable[str]) -> dict[str, str]:
    """Return installed versions for Python and the requested packages."""
    versions = {"python": sys.version.split()[0]}
    for package_name in sorted(set(package_names)):
        try:
            versions[package_name] = importlib.metadata.version(package_name)
        except importlib.metadata.PackageNotFoundError:
            continue
    return versions


# from https://slurm.schedmd.com/job_array.html#env_vars, lower-cased sans SLURM_ prefix
SLURM_KEYS: Final[tuple[str, ...]] = (
    "job_id",
    "array_job_id",
    "array_task_id",
    "array_task_count",
    "mem_per_node",
    "job_nodelist",
    "submit_host",
    "job_partition",
    "job_user",
    "job_account",
    "tasks_per_node",
    "job_qos",
)
SLURM_SUBMIT_KEY: Final[str] = "slurm-submit"


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


def material_id_hash(material_ids: Sequence[str]) -> str:
    """Return an order-independent digest for a material-ID collection."""
    joined_ids = "\n".join(sorted(map(str, material_ids)))
    return hashlib.sha256(joined_ids.encode()).hexdigest()


def material_id_from_atoms(atoms: Atoms) -> str:
    """Extract the material ID from an ASE Atoms info mapping."""
    material_id = atoms.info.get(str(Key.mat_id))
    if not isinstance(material_id, str) or not material_id:
        raise ValueError(f"No {Key.mat_id} in Atoms.info keys {sorted(atoms.info)}")
    return material_id


def _get_calling_file_path(frame: int = 1) -> str:
    """Return the absolute path of the file ``frame`` levels up the call stack."""
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
        account (str): Account to charge for this job.
        partition (str, optional): Slurm partition.
        py_file_path (str, optional): Path to the python script to be submitted.
            Defaults to the path of the file calling slurm_submit().
        slurm_flags (str | list[str], optional): Extra slurm CLI flags. Defaults to ().
            Examples: ('--nodes 1', '--gpus-per-node 1') or ('--mem', '16G').
        array (str, optional): Slurm array specifier. Defaults to None. Example:
            '0-9', '1-10' or '1-10%2'.
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

    if submit_as_temp_file and SLURM_SUBMIT_KEY in sys.argv:
        temp_dir = tempfile.mkdtemp(prefix="slurm_job_")
        temp_file_path = f"{temp_dir}/{os.path.basename(py_file_path)}"
        shutil.copy2(py_file_path, temp_file_path)
        py_file_path = temp_file_path

    if pre_cmd and not pre_cmd.strip().endswith(";"):
        pre_cmd += ";"

    cmd = [
        *("sbatch", "--job-name", job_name),
        *("--output", f"{out_dir}/slurm-%A{'-%a' if array else ''}.log"),
        *(slurm_flags.split() if isinstance(slurm_flags, str) else slurm_flags),
        *("--wrap", f"{pre_cmd} python {py_file_path}".strip()),
    ]
    optional_flags = dict(time=time, account=account, partition=partition, array=array)
    for key, val in optional_flags.items():
        if val is not None:
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
    if pre_cmd:
        slurm_vars["pre_cmd"] = pre_cmd

    # print sbatch cmd into slurm log and at submission time, not to terminal/Jupyter
    if (is_slurm_job and is_log_file) or SLURM_SUBMIT_KEY in sys.argv:
        print(f"\n{' '.join(cmd)}\n".replace(" --", "\n  --"))
    if is_slurm_job and is_log_file:
        for key, val in slurm_vars.items():
            print(f"{key}={val}")

    if SLURM_SUBMIT_KEY not in sys.argv:
        return slurm_vars  # if not submitting slurm job, resume outside code as normal

    result = subprocess.run(cmd, check=True)

    raise SystemExit(result.returncode)


def _is_finite_number(value: object) -> TypeIs[int | float]:
    """Whether value is a usable numeric cost field (NaN would make sums NaN and
    max() depend on unit order).
    """
    return (
        isinstance(value, int | float)
        and not isinstance(value, bool)
        and math.isfinite(value)
    )


def _uniform_label(values: Sequence[object]) -> object | None:
    """The one value every unit shares (compared as sorted JSON), else None.

    Missing (None/empty-string) entries count as disagreement: a label known for
    only some units cannot describe the merged run.
    """
    if not values or any(value is None or value == "" for value in values):
        return None
    normalized = {json.dumps(value, sort_keys=True, default=str) for value in values}
    return values[0] if len(normalized) == 1 else None


def merge_audit_metadata(
    units: Sequence[Mapping[str, object]], *, strict: bool = False
) -> dict[str, object]:
    """Merge per-unit cost and audit provenance into one run's metadata.

    Units (shards, rollout sessions, element files, materials) ran in parallel or in
    sequence, so run/load times and process counts are summed (~ one serial sweep;
    a process charges its load to one unit, see LOAD_CHARGED), memory peaks maxed.
    Cost is all-or-nothing and never blocks merging predictions: a numeric field
    missing from any unit is dropped, a mixed or partial run label drops the whole
    cost block (with a warning), and GPU-dependent costs need a hardware label.

    Raises:
        ValueError: If strict and the units' package ``versions`` differ or are
            missing (not one reproducible prediction set).
    """
    merged: dict[str, object] = {}
    versions = [unit.get("versions") for unit in units]
    if all(isinstance(value, dict) for value in versions) and (
        uniform_versions := _uniform_label(versions)
    ):
        merged["versions"] = uniform_versions
    elif strict:
        raise ValueError("Run units contain mixed or missing package versions")

    dropped: list[str] = []
    for label_key in ("cost_source", *RUN_LABEL_KEYS, "n_shards"):
        values = [unit.get(label_key) for unit in units]
        if (label := _uniform_label(values)) is not None:
            merged[label_key] = label
        elif any(value is not None for value in values):
            distinct = sorted({json.dumps(value, default=str) for value in values})
            dropped.append(f"mixed or partial {label_key} {distinct}")
    if labels_mixed := bool(dropped):
        merged = {key: merged[key] for key in ("versions",) if key in merged}
    numeric_keys = ("run_time_sec", "load_time_sec", "n_processes")
    for key in () if labels_mixed else (*numeric_keys, "max_rss_gb", "max_gpu_mem_gb"):
        values = [value for unit in units if _is_finite_number(value := unit.get(key))]
        if len(values) not in (0, len(units)):
            dropped.append(f"partial {key} ({len(values)}/{len(units)} units)")
        elif values and "hardware" not in merged and key in GPU_DEPENDENT_COST_KEYS:
            dropped.append(f"{key} without a hardware label")
        elif values:
            merged[key] = round(sum(values), 2) if key in numeric_keys else max(values)
    if dropped:
        scope = "the whole cost block" if labels_mixed else "affected cost fields"
        warnings.warn(
            f"Dropping {scope} of {len(units)} merged units: {'; '.join(dropped)}",
            stacklevel=2,
        )

    if completed_at := [
        str(unit["completed_at"]) for unit in units if unit.get("completed_at")
    ]:
        merged["completed_at"] = max(completed_at)
    for source_key in ("hostname", "slurm_job_id", "slurm_array_task_id"):
        if values := sorted(
            {str(unit[source_key]) for unit in units if unit.get(source_key)}
        ):
            merged[f"{source_key}s"] = values
    return merged


def cost_block(metadata: Mapping[str, object]) -> dict[str, object]:
    """The YAML cost mapping (COST_KEYS, canonical order, None values omitted) of
    flat run metadata, whose ``cost_source`` becomes ``source``. Legacy metadata
    without a cost_source keeps only its trustworthy LEGACY_ARTIFACT_COST_KEYS.
    """
    flat = {**metadata, "source": metadata.get("cost_source")}
    keys = LEGACY_ARTIFACT_COST_KEYS if flat["source"] is None else COST_KEYS
    return {key: flat[key] for key in keys if flat.get(key) is not None}


# a cost names hardware and run time (merges dropping either leave labels/RSS only)
PUBLISHABLE_COST_KEYS: Final = {"hardware", "run_time_sec"}


def with_cost_block(
    section: Mapping[str, Any],
    artifact_cost: Mapping[str, object] | None,
    *,
    after: str,
    same_artifact: bool,
    context: str,
    existing: Mapping[str, Any] | None = None,
    replace_author_reported: bool = False,
) -> dict[str, Any]:
    """Return ``section`` with its ``cost`` replaced, as one flow-style YAML line
    right after key ``after`` (first when absent), keeping every other key's order.

    Cost is never mixed from two sources: it is derived wholly from the artifact's
    flat run metadata ``artifact_cost``, or the YAML-only cost of ``existing`` (the
    section before this write, default ``section``) is kept, or there is none.

    1. The artifact records publishable cost: it wins. Recomputing the same artifact
       (``same_artifact``) keeps the YAML's curated ``source`` (an author-run
       artifact stays author_reported); registering a new artifact over an
       author_reported cost needs ``replace_author_reported``.
    2. Otherwise a YAML-only cost (author_reported/salvaged) is kept when
       recomputing the same artifact and dropped (with a warning) for a new one; a
       runner cost is dropped for a new artifact but an error on recompute.

    Raises:
        ValueError: On an undeterminable or invalid source, YAML cost the artifact
            does not back, or an unconfirmed author_reported replacement (messages
            start with ``context``, e.g. ``"mace_mp_0 metrics.md"``).
    """
    cost = _resolve_cost(
        section if existing is None else existing,
        artifact_cost,
        same_artifact=same_artifact,
        context=context,
        replace_author_reported=replace_author_reported,
    )
    out: dict[str, Any] = {} if after in section else dict(cost)
    for key, value in section.items():
        if key != "cost":
            out[key] = value
        if key == after:
            out |= cost
    return out


def _resolve_cost(
    existing: Mapping[str, Any],
    artifact_cost: Mapping[str, object] | None,
    *,
    same_artifact: bool,
    context: str,
    replace_author_reported: bool,
) -> dict[str, object]:
    """with_cost_block's ``{"cost": <flow mapping>}`` entry, or ``{}`` for none."""
    from ruamel.yaml.comments import CommentedMap

    existing_block = dict(existing.get("cost") or {})
    existing_source = existing_block.get("source")
    new_block = cost_block(artifact_cost or {})
    if existing_block and existing_source not in COST_SOURCES:
        raise ValueError(
            f"{context}: existing cost {existing_block} has no valid source "
            f"(expected one of {COST_SOURCES})"
        )
    if existing_source == "author_reported" and not (
        same_artifact or replace_author_reported
    ):
        raise ValueError(
            f"{context}: refusing to replace author_reported cost {existing_block} "
            f"with artifact cost {new_block}; pass replace_author_reported=True to "
            "confirm"
        )
    if not new_block.keys() >= PUBLISHABLE_COST_KEYS:
        if existing_source in YAML_ONLY_COST_SOURCES and same_artifact:
            return {"cost": existing.get("cost")}
        if existing_block and same_artifact:
            raise ValueError(
                f"{context}: YAML cost {existing_block} is tagged {existing_source!r} "
                "but the registered artifact records no cost"
            )
        if existing_source in YAML_ONLY_COST_SOURCES:
            warnings.warn(
                f"{context}: dropping {existing_source} cost {existing_block}: it "
                "describes the replaced artifact, and the new one records no "
                "publishable cost",
                stacklevel=2,
            )
        return {}
    if same_artifact and existing_source == "author_reported":  # author-run artifact
        new_block["source"] = existing_source
    elif "source" not in new_block and existing_source not in ("runner", "salvaged"):
        raise ValueError(
            f"{context}: artifact cost {new_block} records no cost_source and the "
            f"YAML has no artifact-backed tag ({existing_source=})"
        )
    new_block.setdefault("source", existing_source)
    if new_block["source"] not in COST_SOURCES:
        raise ValueError(
            f"{context}: artifact cost_source {new_block['source']!r} must be one of "
            f"{COST_SOURCES}"
        )
    cost = CommentedMap({key: new_block[key] for key in COST_KEYS if key in new_block})
    cost.fa.set_flow_style()  # one YAML line, units implied by the key suffixes
    return {"cost": cost}


# YAML task sections (dotted under metrics) holding a cost mapping, and the file ref
# whose artifact records that cost
COST_SECTIONS: Final = {
    "discovery": "run_info_file",
    "diatomics": "pred_file",
    "md": "pred_file",
    "phonons.kappa_103": "run_info_file",
    "phonons.kappa_103.harmonic": "run_info_file",
}


def cost_provenance_problems(
    metrics: Mapping[str, object], artifact_root: str
) -> list[str]:
    """Cross-check every COST_SECTIONS ``cost`` in a model's ``metrics`` against its
    artifact's recorded cost (aggregated MD CSV, else ``run_metadata``) when that
    exists under ``artifact_root`` (never downloaded). Publishable artifact cost
    must match key for key (numbers to 1e-6; an author_reported source is curated),
    legacy artifacts (no cost_source) on LEGACY_ARTIFACT_COST_KEYS only. Otherwise
    only YAML-only sources are valid, and tagged cost the YAML omits is flagged.
    """
    import gzip

    from matbench_discovery.data import file_ref_name

    problems: list[str] = []
    for section_path, ref_key in COST_SECTIONS.items():
        section: Any = metrics
        for key in section_path.split("."):
            section = section.get(key) if isinstance(section, Mapping) else None
        if not isinstance(section, Mapping):
            continue
        yaml_cost = dict(section.get("cost") or {})
        yaml_only = yaml_cost.get("source") in YAML_ONLY_COST_SOURCES
        ref_name = file_ref_name(section.get(ref_key))
        artifact_path = f"{artifact_root}/{ref_name}"
        if ref_name is None or not os.path.isfile(artifact_path):
            if ref_name is None and yaml_cost and not yaml_only:
                problems.append(f"{section_path}: {yaml_cost} without a {ref_key}")
            continue
        if section_path == "md":
            import pandas as pd

            from matbench_discovery.metrics.md import aggregate_md_cost

            df_md = pd.read_csv(artifact_path, index_col="system")
            artifact_cost = cost_block(aggregate_md_cost(df_md))
        else:
            open_fn = gzip.open if artifact_path.endswith(".gz") else open
            with open_fn(artifact_path, mode="rt", encoding="utf-8") as file:
                data = json.load(file)
            # diatomics files predating run_metadata record no cost; run-info
            # sidecars exist only to hold it
            default = {} if section_path == "diatomics" else None
            if not isinstance(run_metadata := data.get("run_metadata", default), dict):
                raise TypeError(f"{artifact_path} has no run_metadata mapping")
            artifact_cost = cost_block(run_metadata)
        if not artifact_cost.keys() >= PUBLISHABLE_COST_KEYS:
            if yaml_cost and not yaml_only:
                problems.append(
                    f"{section_path}: cost {yaml_cost} but the artifact records none"
                )
            continue
        if not yaml_cost:
            if "source" in artifact_cost:
                problems.append(f"{section_path}: unregistered cost {artifact_cost}")
            continue
        legacy = "source" not in artifact_cost
        for key in LEGACY_ARTIFACT_COST_KEYS if legacy else COST_KEYS:
            if key == "source" and yaml_cost["source"] == "author_reported":
                continue  # author-run artifact: the YAML tag is curated
            yaml_val, artifact_val = yaml_cost.get(key), artifact_cost.get(key)
            if _is_finite_number(yaml_val) and _is_finite_number(artifact_val):
                if math.isclose(yaml_val, artifact_val, rel_tol=0, abs_tol=1e-6):
                    continue
            elif yaml_val == artifact_val:
                continue
            problems.append(
                f"{section_path}: {key}: YAML {yaml_val!r} != artifact {artifact_val!r}"
            )
    return problems
