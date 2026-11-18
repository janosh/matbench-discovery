"""Molecular dynamics rollout helpers for the MD benchmark task.

Runs NVT simulations with MLIP calculators using the DynaMat v1.0 protocol:
Nose-Hoover chain thermostat, 0.25 fs timestep, 25 fs thermostat time scale, and
frames recorded every 10 steps for 20 ps.
"""

import contextlib
import dataclasses
import functools
import json
import os
import time
import zipfile
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from typing import Any, cast

import h5py
import numpy as np
import pandas as pd
from matscikit import NoseHooverChainNVT, Structure
from matscikit.io import Frame, read_extxyz, write_extxyz
from matscikit.md import maxwell_boltzmann_momenta
from matscikit.model import Model as ForceModel
from tqdm import tqdm

from matbench_discovery import (
    MATSCIKIT_EXTENSION,
    atomic_output_path,
    file_digest,
    today,
)
from matbench_discovery.hpc import (
    LOAD_CHARGED,
    check_resume_labels,
    memory_peaks,
    merge_audit_metadata,
    reset_gpu_peak_memory,
)
from matbench_discovery.metrics import md as md_metrics
from matbench_discovery.trajectory import Trajectory

# DynaMat v1.0 protocol thermostat time scale (fs) shared by every rollout entry point
THERMOSTAT_TIME_SCALE_FS = 25.0

# === single-file multi-system reference dataset ===
# DynaMat v1.0 references ship as one HDF5 with a group per system, each carrying the
# trajectory arrays (guarded by Trajectory's per-group TRAJECTORY_SCHEMA) plus its
# saved-frame interval (dt_fs) and target temperature as group attrs. Self-describing
# -> no separate settings CSV or system-name matching.


def default_md_reference_path() -> str:
    """Path to the auto-downloaded DynaMat v1.0 reference HDF5."""
    from matbench_discovery.enums import DataFiles

    return DataFiles.dynamat_v1_0_md_trajectories.path


def write_reference_h5(
    path: str, entries: Iterable[tuple[str, Trajectory, float, float]]
) -> None:
    """Atomically write a multi-system reference HDF5: one group per system holding the
    label-free trajectory arrays plus ``dt_fs`` (saved-frame interval) and
    ``temperature_kelvin`` attrs. ``entries`` yields (system_name, trajectory, dt_fs,
    temperature_kelvin); it may be a lazy generator so the converter parses one
    trajectory at a time and writes it straight to disk, bounding peak memory to the
    single largest system.
    """
    n_written = 0
    with atomic_output_path(path) as tmp_path, h5py.File(tmp_path, "w") as file:
        for system_name, trajectory, dt_fs, temperature_kelvin in entries:
            group = file.create_group(system_name)
            # strip supervised labels so the public artifact can't double as a
            # validation set; dataclasses.replace keeps any future fields intact
            public_trajectory = dataclasses.replace(
                trajectory, energy=None, forces=None
            )
            public_trajectory.write_to_h5_group(group)
            group.attrs["dt_fs"] = float(dt_fs)
            group.attrs["temperature_kelvin"] = float(temperature_kelvin)
            n_written += 1
        if n_written == 0:  # nothing yielded: the empty temp file is dropped
            raise ValueError("Cannot write a reference file with no systems")


def list_reference_systems(path: str) -> list[str]:
    """Sorted system names (group keys) in a reference HDF5 file."""
    with h5py.File(path, "r") as file:
        return sorted(name for name in file if isinstance(file[name], h5py.Group))


def read_reference_trajectory(
    path: str, system_name: str, *, max_frames: int | None = None
) -> tuple[Trajectory, float, float]:
    """Read one system from a reference HDF5, returning (trajectory, dt_fs,
    temperature_kelvin). Only the requested frames are read/decompressed, so dry-run
    head reads (``max_frames``) of huge references stay fast.
    """
    frames = slice(0, max_frames) if max_frames is not None else slice(None)
    with h5py.File(path, "r") as file:
        if system_name not in file:
            raise KeyError(
                f"{system_name!r} not in reference {path!r}; e.g. {sorted(file)[:5]}"
            )
        group = file[system_name]
        trajectory = Trajectory.read_from_h5_group(group, frames=frames)
        dt_fs = float(group.attrs["dt_fs"])
        temperature_kelvin = float(group.attrs["temperature_kelvin"])
    return trajectory, dt_fs, temperature_kelvin


def read_trajectory(file_path: str, *, index: int | slice = slice(None)) -> list[Frame]:
    """Frames of an (optionally compressed) extxyz trajectory file; ``index`` selects
    frames like a list index (an int returns a one-frame list).
    """
    return read_extxyz(file_path, index=index)


def run_nvt_md(
    atoms: Structure,
    model: ForceModel,
    *,
    temperature_kelvin: float,
    n_steps: int = 80_000,
    time_step_fs: float = 0.25,
    record_interval: int = 10,
    thermostat_time_scale_fs: float = THERMOSTAT_TIME_SCALE_FS,
    seed: int = 0,
    progress_interval: int = 1_000,
    progress_label: str = "MD rollout",
) -> list[Frame]:
    """Run one NVT MD simulation and return recorded frames.

    Velocities are initialized from a Maxwell-Boltzmann distribution. Each recorded
    frame carries energy, forces and (when the model supports it) stress plus momenta,
    so trajectories written to extxyz retain everything needed for RDF/vDOS/pressure
    evaluation.

    Args:
        atoms: Initial structure (not modified).
        model: ForceModel providing energies, forces and ideally stress.
        temperature_kelvin: Target temperature in Kelvin.
        n_steps: Number of MD steps. Default 80,000 = 20 ps at 0.25 fs.
        time_step_fs: Integration time step in femtoseconds.
        record_interval: Record a frame every this many steps.
        thermostat_time_scale_fs: Nose-Hoover damping time scale in femtoseconds.
        seed: Seed for the Maxwell-Boltzmann velocity initialization.
        progress_interval: Update the progress bar every this many steps (useful for
            monitoring long rollouts in batch job logs). 0 disables the progress bar.
        progress_label: Description shown beside the progress bar.

    Returns:
        list[Frame]: Recorded frames including the initial one, i.e.
            n_steps // record_interval + 1 frames.
    """
    dynamics = _thermalized_nvt_dynamics(
        atoms,
        model,
        temperature_kelvin=temperature_kelvin,
        time_step_fs=time_step_fs,
        thermostat_time_scale_fs=thermostat_time_scale_fs,
        seed=seed,
    )
    frames: list[Frame] = [_snapshot_frame(dynamics.frame(0))]
    _run_dynamics_with_progress(
        dynamics,
        model,
        n_steps=n_steps,
        record_interval=record_interval,
        on_frame=lambda frame: frames.append(_snapshot_frame(frame)),
        interval=progress_interval,
        label=progress_label,
    )
    return frames


def _supports_stress(model: ForceModel, structure: Structure) -> bool:
    """Whether ``model`` returns stress (probed once on the initial structure).

    Only NotImplementedError (ASE's PropertyNotImplementedError subclasses it) means
    "no stress head"; any other failure is a real model error and propagates.
    """
    try:
        (result,) = model.compute([structure], ["energy", "forces", "stress"])
    except NotImplementedError:
        return False
    return result.get("stress") is not None


def _thermalized_nvt_dynamics(
    atoms: Structure,
    model: ForceModel,
    *,
    temperature_kelvin: float,
    time_step_fs: float,
    thermostat_time_scale_fs: float,
    seed: int,
) -> NoseHooverChainNVT:
    """Seed Maxwell-Boltzmann momenta and build the Nose-Hoover chain integrator
    (evaluated at step 0). Single deterministic setup shared by the in-memory and
    resumable rollouts so both thermalize identically for a given seed.
    """
    momenta = maxwell_boltzmann_momenta(
        atoms.masses, temperature_kelvin, np.random.default_rng(seed)
    )
    dynamics = NoseHooverChainNVT(
        [atoms],
        [momenta],
        temperature_k=temperature_kelvin,
        time_step_fs=time_step_fs,
        damping_fs=thermostat_time_scale_fs,
        stress=_supports_stress(model, atoms),
    )
    dynamics.run(model, 0)  # evaluate the initial state
    return dynamics


def _snapshot_frame(state: dict[str, Any]) -> Frame:
    """Frame of an integrator state dict with energy/forces/(stress) and momenta."""
    structure = state["structure"].copy()
    structure.info["md_step"] = int(state["step"])
    return Frame(
        structure,
        arrays={"momenta": np.asarray(state["momenta"])},
        energy=float(state["energy"]),
        forces=np.asarray(state["forces"]),
        stress=None if state.get("stress") is None else np.asarray(state["stress"]),
    )


def _run_dynamics_with_progress(
    dynamics: NoseHooverChainNVT,
    model: ForceModel,
    *,
    n_steps: int,
    record_interval: int,
    on_frame: Callable[[dict[str, Any]], None],
    interval: int,
    label: str,
) -> None:
    """Advance to ``n_steps`` total steps, calling ``on_frame`` every
    ``record_interval`` steps, with a tqdm progress bar (0 interval disables it).
    """
    done = int(dynamics.get_state(0)["n_steps"])
    remaining_steps = n_steps - done
    if remaining_steps <= 0:
        return

    def callback(_idx: int, state: dict[str, Any]) -> None:
        on_frame(state)

    if interval <= 0:
        dynamics.run(
            model, remaining_steps, callback=callback, interval=record_interval
        )
        return
    with tqdm(
        total=n_steps, initial=done, desc=label, unit="step", mininterval=5, leave=True
    ) as pbar:
        while remaining_steps > 0:
            chunk = min(interval, remaining_steps)
            dynamics.run(model, chunk, callback=callback, interval=record_interval)
            pbar.update(chunk)
            remaining_steps -= chunk


@dataclass(frozen=True)
class NvtParams:
    """The rollout parameters a checkpoint must match to be safely resumable, i.e.
    those defining the deterministic trajectory and recording cadence. ``n_steps`` (run
    length) is deliberately excluded: a checkpoint at step N resumes for any target
    length >= N, so it's bound-checked rather than matched for equality.
    """

    temperature_kelvin: float
    time_step_fs: float
    record_interval: int
    thermostat_time_scale_fs: float
    seed: int


# === crash-safe checkpointing for resumable rollouts ===
# NoseHooverChainNVT is deterministic after the seeded velocity init (no per-step RNG),
# so its full state is positions, momenta, thermostat (eta, p_eta) and the step count;
# restoring them reproduces an uninterrupted run to ~1e-14 (positions pass through
# fractional storage, ~1 ulp, which grows only mildly over a rollout). Schema-tagged and
# matscikit-version-gated so an integrator change fails closed (recompute) rather than
# silently corrupting trajectories.
# Schema 3: the cost units of the sessions so far (runtime + steady-state run time +
# memory peaks), whose run labels a resume must match (never mixing GPUs); schema 4
# moved the integrator from ASE to matscikit.
CHECKPOINT_SCHEMA = 4


@functools.cache
def matscikit_build_digest() -> str:
    """SHA-256 of the compiled matscikit extension, so a checkpoint written by other
    numerics (``matscikit.__version__`` is unchanged across rebuilds) is never resumed.
    """
    return file_digest(MATSCIKIT_EXTENSION)


NpzValue = np.ndarray | np.generic | str | int | float | bool


def run_info_path(traj_path: str) -> str:
    """Path of the JSON sidecar storing rollout cost and hardware provenance.

    Written next to the trajectory once a rollout completes so later jobs that reuse
    the cached rollout (e.g. metric recomputes) can still report its cost. Trajectory
    files carry no timestamps, so this data is unrecoverable once lost - hence a
    persistent sidecar rather than an in-memory return value.
    """
    return f"{traj_path.removesuffix('.extxyz')}-run-info.json"


def read_run_info(traj_path: str) -> dict[str, object]:
    """Load a rollout's run-info sidecar, or an empty dict when absent (rollouts that
    predate sidecars). A corrupt sidecar raises: silently dropping it would publish
    the rollout without its cost provenance or reuse it without a parameter check.
    """
    path = run_info_path(traj_path)
    if not os.path.isfile(path):
        return {}
    with open(path, encoding="utf-8") as file:
        try:
            run_info = json.load(file)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Corrupt MD run-info sidecar {path}: {exc}") from exc
    if not isinstance(run_info, dict):
        raise TypeError(f"MD run-info sidecar {path} must hold a JSON object")
    return run_info


def run_info_row(run_info: Mapping[str, object]) -> dict[str, float | int | str]:
    """A run-info sidecar's cost fields as per-system CSV columns. ``acceleration``
    is ``+``-joined, ``none`` for [] (decoded by metrics.md.aggregate_md_cost).
    """
    row: dict[str, float | int | str] = {}
    for key in md_metrics.MD_COST_COLS:
        value = run_info.get(key)
        if key == "acceleration" and isinstance(value, list):
            row[key] = "+".join(map(str, value)) or "none"
        elif isinstance(value, str | int | float) and not isinstance(value, bool):
            row[key] = value
    return row


def _atomic_savez(path: str, **arrays: NpzValue) -> None:
    """Write an NPZ file through a temp path, then atomically promote it."""
    with atomic_output_path(path) as tmp_path, open(tmp_path, "wb") as file:
        np.savez(file, **cast("dict[str, Any]", arrays))


def _nvt_state(state: dict[str, Any]) -> dict[str, np.ndarray]:
    """Full integrator state (positions, momenta, eta, p_eta, n_steps) of a frame
    dict, as passed to ``NoseHooverChainNVT.run`` callbacks or from ``frame()``.
    """
    return {
        "positions": np.asarray(state["structure"].positions),
        "momenta": np.asarray(state["momenta"]),
        "eta": np.asarray(state["eta"]),
        "p_eta": np.asarray(state["p_eta"]),
        "n_steps": np.asarray(int(state["step"])),
    }


def _load_checkpoint(
    ckpt_path: str, *, atoms: Structure, params: NvtParams, n_steps: int
) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]] | None:
    """Load a checkpoint with its sessions' cost units if present,
    readable and schema/matscikit-build/parameter-matched; else None (so the caller
    recomputes from scratch rather than resuming a stale/foreign state).
    """
    if not os.path.isfile(ckpt_path):
        return None
    try:
        with np.load(ckpt_path, allow_pickle=False) as data:
            ckpt = {key: data[key] for key in data.files}
    except (OSError, ValueError, EOFError, zipfile.BadZipFile) as exc:
        print(f"  ignoring unreadable checkpoint {ckpt_path}: {exc}")
        return None

    target_frames = n_steps // params.record_interval + 1
    try:
        schema = int(ckpt["schema"])
        engine_build = str(ckpt["matscikit_build"])
    except KeyError, TypeError, ValueError:
        print(f"  ignoring checkpoint {ckpt_path}: missing/corrupt schema metadata")
        return None
    if schema != CHECKPOINT_SCHEMA or engine_build != matscikit_build_digest():
        print(f"  ignoring checkpoint {ckpt_path}: schema/matscikit build mismatch")
        return None

    try:
        sessions = json.loads(str(ckpt["sessions"]))
        stored = NvtParams(
            temperature_kelvin=float(ckpt["temperature_kelvin"]),
            time_step_fs=float(ckpt["time_step_fs"]),
            record_interval=int(ckpt["record_interval"]),
            thermostat_time_scale_fs=float(ckpt["thermostat_time_scale_fs"]),
            seed=int(ckpt["seed"]),
        )
        nsteps = int(ckpt["nsteps"])
        n_frames = int(ckpt["n_frames"])
        consistent = (
            stored == params
            and nsteps <= n_steps
            and n_frames <= target_frames
            and isinstance(sessions, list)
            and all(isinstance(session, dict) for session in sessions)
            and nsteps == (n_frames - 1) * params.record_interval
            and list(ckpt["symbols"]) == atoms.symbols
            and np.allclose(ckpt["cell"], atoms.lattice)
            and np.array_equal(ckpt["pbc"], atoms.pbc)
        )
    except KeyError, TypeError, ValueError:  # JSONDecodeError is a ValueError
        consistent = False
    if not consistent:
        print(f"  ignoring checkpoint {ckpt_path}: schema/version/parameters mismatch")
        return None
    return ckpt, sessions


def _restore_nvt_state(
    dynamics: NoseHooverChainNVT, ckpt: dict[str, np.ndarray]
) -> None:
    """Restore integrator + thermostat state from a loaded checkpoint."""
    dynamics.set_state(
        0,
        ckpt["q"],
        ckpt["p"],
        int(ckpt["nsteps"]),
        eta=ckpt["eta"].tolist(),
        p_eta=ckpt["p_eta"].tolist(),
    )


def run_nvt_md_resumable(
    atoms: Structure,
    model: ForceModel,
    *,
    out_path: str,
    temperature_kelvin: float,
    n_steps: int = 80_000,
    time_step_fs: float = 0.25,
    record_interval: int = 10,
    thermostat_time_scale_fs: float = THERMOSTAT_TIME_SCALE_FS,
    seed: int = 0,
    progress_interval: int = 1_000,
    progress_label: str = "MD rollout",
    checkpoint_every_n_frames: int = 200,
    runtime: Mapping[str, Any],
) -> list[Frame]:
    """Run an NVT rollout that survives interruption (timeouts/preemption).

    Frames stream to ``{out_path}.part`` and the integrator state is checkpointed to
    ``{out_path}.ckpt.npz`` every ``checkpoint_every_n_frames`` frames. A rerun resumes
    from the last checkpoint (truncating any un-checkpointed tail), continuing the same
    deterministic trajectory. Only on completion is ``.part`` atomically promoted to
    ``out_path`` and the checkpoint removed. Returns the full trajectory.

    Each session's cost unit (``runtime`` from calculators.load_calculator_with_runtime,
    steady-state run time, memory peaks) persists in the checkpoint and, on
    completion, all are merged into a ``run_info_path(out_path)`` sidecar, so the
    rollout's cost survives later reuse.

    Raises:
        ValueError: If a checkpoint written under other run labels would be resumed.
    """
    # .part.extxyz (not .part) so the partial file is still recognizably extxyz
    part_path = f"{out_path.removesuffix('.extxyz')}.part.extxyz"
    ckpt_path = f"{out_path}.ckpt.npz"
    target_frames = n_steps // record_interval + 1
    if checkpoint_every_n_frames < 1:
        raise ValueError(
            f"{checkpoint_every_n_frames=} must be positive for MD checkpointing"
        )

    dynamics = _thermalized_nvt_dynamics(
        atoms,
        model,
        temperature_kelvin=temperature_kelvin,
        time_step_fs=time_step_fs,
        thermostat_time_scale_fs=thermostat_time_scale_fs,
        seed=seed,
    )
    params = NvtParams(
        temperature_kelvin=temperature_kelvin,
        time_step_fs=time_step_fs,
        record_interval=record_interval,
        thermostat_time_scale_fs=thermostat_time_scale_fs,
        seed=seed,
    )

    loaded = _load_checkpoint(ckpt_path, atoms=atoms, params=params, n_steps=n_steps)
    frames_written = 0
    sessions: list[dict[str, Any]] = []  # cost units of earlier (timed-out) sessions
    # attribute peak VRAM to this rollout even when one process runs several systems
    reset_gpu_peak_memory(runtime)
    if loaded is not None and os.path.isfile(part_path):
        ckpt, stored_sessions = loaded
        n_ckpt = int(ckpt["n_frames"])
        # keep only checkpointed frames; a kill may have appended an un-checkpointed
        # (possibly malformed) tail, so read just the first n_ckpt and rewrite cleanly
        try:
            kept = read_trajectory(part_path, index=slice(None, n_ckpt))
        except (OSError, ValueError, IndexError) as exc:
            print(f"  resume failed reading {part_path} ({exc}); recomputing")
            kept = []
        if len(kept) == n_ckpt:
            check_resume_labels(stored_sessions[-1], runtime, ckpt_path)
            sessions = stored_sessions
            _restore_nvt_state(dynamics, ckpt)
            dynamics.run(model, 0)  # re-evaluate the restored state
            write_extxyz(part_path, kept)  # truncate to the checkpointed frames
            frames_written = n_ckpt
            print(f"  resuming {out_path} from frame {frames_written}/{target_frames}")
    session_start = time.perf_counter()

    def session_unit() -> dict[str, Any]:
        """This session's cost unit so far (after MD steps ran)."""
        run_time_sec = time.perf_counter() - session_start
        return {**runtime, "run_time_sec": run_time_sec} | memory_peaks(
            runtime, succeeded=True
        )

    if frames_written == 0:  # fresh start: clear any stale partial files
        for stale in (part_path, ckpt_path):
            with contextlib.suppress(FileNotFoundError):
                os.remove(stale)

    def record_and_checkpoint(state: dict[str, Any]) -> None:
        """Append the current frame and periodically checkpoint the NVT state."""
        nonlocal frames_written
        write_extxyz(part_path, _snapshot_frame(state), append=True)
        frames_written += 1
        if frames_written % checkpoint_every_n_frames == 0:
            nvt_state = _nvt_state(state)
            _atomic_savez(
                ckpt_path,
                schema=CHECKPOINT_SCHEMA,
                matscikit_build=matscikit_build_digest(),
                nsteps=int(nvt_state["n_steps"]),
                n_steps=n_steps,
                n_frames=frames_written,
                q=nvt_state["positions"],
                p=nvt_state["momenta"],
                eta=nvt_state["eta"],
                p_eta=nvt_state["p_eta"],
                record_interval=params.record_interval,
                time_step_fs=params.time_step_fs,
                temperature_kelvin=params.temperature_kelvin,
                thermostat_time_scale_fs=params.thermostat_time_scale_fs,
                seed=params.seed,
                symbols=np.array(atoms.symbols),
                cell=atoms.lattice,
                pbc=np.array(atoms.pbc),
                # so a killed session's cost still counts after resume
                sessions=json.dumps([*sessions, session_unit()]),
            )

    if frames_written == 0:  # record the initial state (step 0)
        record_and_checkpoint(dynamics.frame(0))
    _run_dynamics_with_progress(
        dynamics,
        model,
        n_steps=n_steps,
        record_interval=record_interval,
        on_frame=record_and_checkpoint,
        interval=progress_interval,
        label=progress_label,
    )

    frames = read_trajectory(part_path)
    if len(frames) != target_frames:
        raise RuntimeError(
            f"{out_path} rollout produced {len(frames)} frames != {target_frames}"
        )
    # cost keys only from the merge, never one session's values; nvt_params make a
    # later run with other settings recompute (validate_pred_trajectory)
    cost_cols = md_metrics.MD_COST_COLS
    run_info = {key: val for key, val in runtime.items() if key not in cost_cols} | {
        **merge_audit_metadata([*sessions, session_unit()], strict=True),
        "nvt_params": dataclasses.asdict(params),
        "completed_at": datetime.now(UTC).isoformat(),
    }
    with (
        atomic_output_path(run_info_path(out_path)) as tmp_path,
        open(tmp_path, mode="w", encoding="utf-8") as file,
    ):
        json.dump(run_info, file, indent=2)
    os.replace(part_path, out_path)  # atomic promote only once complete
    with contextlib.suppress(FileNotFoundError):
        os.remove(ckpt_path)
    return frames


def validate_pred_trajectory(
    path: str, *, expected_frames: int, ref_atoms: Structure, params: NvtParams
) -> list[Frame] | None:
    """Read and validate a cached predicted rollout, returning its frames if sound or
    None if it should be recomputed.

    A truncated or otherwise inconsistent rollout would silently evaluate as a shorter
    trajectory (evaluate_md_system time-matches to the shorter span), so reuse only
    when the file is complete and consistent: run with ``params`` (when its run-info
    sidecar records them; rollouts predating that record are trusted), readable,
    exactly ``expected_frames`` frames with a strict ``md_step`` sequence 0,
    record_interval, 2*record_interval, ...; matching atom count and chemical symbols,
    cell and pbc of ``ref_atoms``; and all-finite positions and velocities.
    """
    recorded_params = read_run_info(path).get("nvt_params")
    if isinstance(recorded_params, dict) and NvtParams(**recorded_params) != params:
        print(f"  recomputing {path}: run with {recorded_params} != {params}")
        return None
    record_interval = params.record_interval
    try:
        frames = read_trajectory(path)
    except (OSError, ValueError, IndexError, StopIteration) as exc:
        print(f"  recomputing {path}: unreadable ({exc})")
        return None

    if (n_frames := len(frames)) != expected_frames:
        print(f"  recomputing {path}: {n_frames} frames != expected {expected_frames}")
        return None

    ref_symbols = ref_atoms.symbols
    for idx, frame in enumerate(frames):
        structure = frame.structure
        momenta = frame.arrays.get("momenta")
        reason = None
        if int(structure.info.get("md_step", -1)) != idx * record_interval:
            reason = (
                f"md_step {structure.info.get('md_step')} != {idx * record_interval}"
            )
        elif structure.symbols != ref_symbols:
            reason = "chemical symbols differ from reference"
        elif (
            not np.allclose(structure.lattice, ref_atoms.lattice)
            or structure.pbc != ref_atoms.pbc
        ):
            reason = "cell or pbc differ from reference"
        elif not np.isfinite(structure.positions).all() or (
            momenta is not None and not np.isfinite(momenta).all()
        ):
            reason = "non-finite positions or velocities"
        if reason:
            print(f"  recomputing {path}: frame {idx} {reason}")
            return None

    return frames


def run_md_benchmark(
    *,
    model: ForceModel,
    model_key: str,
    out_dir: str,
    ref_file: str | None = None,
    private_ref_file: str | None = None,
    systems: list[str] | None = None,
    n_steps: int = 80_000,
    time_step_fs: float = 0.25,
    record_interval: int = 10,
    seed: int = 0,
    dry_run: bool = False,
    runtime: Mapping[str, Any],
) -> pd.DataFrame:
    """Run NVT rollouts and compute MD metrics for one model across reference systems.

    For each system in the reference HDF5, this rolls out an NVT trajectory from the
    reference initial structure (reusing an existing rollout if present), then compares
    it to the reference via RDF, ADF, vDOS and pressure observables. If a private
    labeled reference file is provided, diagnostic energy/force RMSEs are computed from
    that file's hidden labels. Writes one gzipped per-system metrics CSV and returns
    the per-system DataFrame. Each row also carries ``n_atoms`` plus the rollout's cost
    provenance (metrics.md.MD_COST_COLS) read from the trajectory's run-info sidecar,
    which survives metric-recompute jobs that reuse cached rollouts.

    Args:
        model: ForceModel under test (matscikit ``compute`` protocol).
        model_key: ForceModel enum name/key, used in output filenames.
        out_dir: Directory for rollout trajectories and the metrics CSV.
        ref_file: Reference HDF5 (one group per system, each carrying dt_fs and
            temperature_kelvin attrs). Defaults to the auto-downloaded DynaMat v1.0 set.
        private_ref_file: Optional private HDF5 with matching system groups and
            energy/force labels for maintainer-only diagnostic RMSEs.
        systems: Subset of system names to run. Defaults to all in ref_file.
        n_steps: MD steps per rollout (80,000 = 20 ps at 0.25 fs).
        time_step_fs: MD integration time step.
        record_interval: Record a frame every this many MD steps.
        seed: Maxwell-Boltzmann velocity seed.
        dry_run: Smoke test: one system, a handful of MD steps, and a capped reference
            slice, to verify the model + pipeline run error-free in seconds without
            writing outputs. Returned metrics are not meaningful.
        runtime: calculators.load_calculator_with_runtime() provenance, recorded in
            each new rollout's run info.

    Returns:
        pd.DataFrame: One row per system indexed by system name.
    """
    if ref_file is None:
        ref_file = default_md_reference_path()

    all_systems = list_reference_systems(ref_file)
    if systems is not None and (unknown := sorted(set(systems) - set(all_systems))):
        raise ValueError(
            f"Unknown systems {unknown} not in reference {ref_file!r}, available: "
            f"{all_systems}"
        )
    system_names = (
        all_systems
        if systems is None
        else [name for name in all_systems if name in systems]
    )
    if not system_names:  # clearer than the downstream empty-DataFrame KeyError
        raise ValueError(f"No systems found in reference {ref_file!r} for {systems=}")
    if dry_run:  # one system, a few steps, capped reference slice
        system_names = system_names[:1]
        # enough recorded frames that time-matching keeps >=4 even for the coarsest
        # reference cadence (dt_fs up to 2 fs vs the 2.5 fs prediction cadence)
        n_steps = record_interval * 20
    else:
        os.makedirs(out_dir, exist_ok=True)  # for rollout trajectories and metrics CSV

    pred_time_step_fs = time_step_fs * record_interval
    rows: list[dict[str, float | int | str]] = []
    for system_name in tqdm(system_names, desc=f"MD systems ({model_key})"):
        progress_label = f"{model_key} {system_name}"
        # dry run reads only the first 64 reference frames; the HDF5 reader only
        # decompresses the frames it returns
        ref_trajectory, ref_time_step_fs, temperature_kelvin = (
            read_reference_trajectory(
                ref_file, system_name, max_frames=64 if dry_run else None
            )
        )
        # abort before any GPU time is spent if the reference's dt_fs attr is stale
        md_metrics.check_reference_dt(
            ref_trajectory[:256],
            time_step_fs=ref_time_step_fs,
            temperature_kelvin=temperature_kelvin,
            label=f"reference {system_name!r} in {ref_file!r}",
        )
        initial_atoms = ref_trajectory.frame_as_structure(0)

        pred_traj_path = f"{out_dir}/{system_name}-nvt-{model_key}.extxyz"
        target_frames = n_steps // record_interval + 1
        if dry_run:  # in-memory run, no checkpoint/output files
            pred_frames = run_nvt_md(
                initial_atoms,
                model,
                temperature_kelvin=temperature_kelvin,
                n_steps=n_steps,
                time_step_fs=time_step_fs,
                record_interval=record_interval,
                seed=seed,
                progress_label=progress_label,
            )
        else:
            # reuse only a complete, consistent rollout; a truncated file would
            # otherwise silently evaluate as a shorter trajectory
            pred_frames = None
            if os.path.isfile(pred_traj_path):
                pred_frames = validate_pred_trajectory(
                    pred_traj_path,
                    expected_frames=target_frames,
                    ref_atoms=initial_atoms,
                    params=NvtParams(
                        temperature_kelvin=temperature_kelvin,
                        time_step_fs=time_step_fs,
                        record_interval=record_interval,
                        thermostat_time_scale_fs=THERMOSTAT_TIME_SCALE_FS,
                        seed=seed,
                    ),
                )
            if pred_frames is None:  # crash-safe resumable rollout to pred_traj_path
                pred_frames = run_nvt_md_resumable(
                    initial_atoms,
                    model,
                    out_path=pred_traj_path,
                    temperature_kelvin=temperature_kelvin,
                    n_steps=n_steps,
                    time_step_fs=time_step_fs,
                    record_interval=record_interval,
                    seed=seed,
                    progress_label=progress_label,
                    runtime=runtime,
                )
                runtime = {**runtime, **LOAD_CHARGED}  # its first rollout carries it

        system_metrics = md_metrics.evaluate_md_system(
            ref_trajectory,
            pred_frames,  # evaluate_md_system coerces the frame list to a Trajectory
            ref_time_step_fs=ref_time_step_fs,
            pred_time_step_fs=pred_time_step_fs,
            progress_label=progress_label,
        )
        if private_ref_file is not None:
            # only the labeled trajectory is needed here; the private file's dt_fs and
            # temperature are already validated on the public reference above
            private_ref_trajectory = read_reference_trajectory(
                private_ref_file, system_name, max_frames=64 if dry_run else None
            )[0]
            # guard against a stale/mismatched private file: matching group names
            # with different structures would yield plausible-but-wrong diagnostics
            if not np.array_equal(
                private_ref_trajectory.atomic_numbers, ref_trajectory.atomic_numbers
            ):
                raise ValueError(
                    f"private reference {system_name!r} has different atoms than the "
                    f"public reference ({private_ref_trajectory.n_atoms} vs "
                    f"{ref_trajectory.n_atoms}); stale {private_ref_file=}?"
                )
            system_metrics |= md_metrics.calc_energy_force_rmse(
                private_ref_trajectory, model
            )
        row: dict[str, float | int | str] = {
            "system": system_name,
            "temperature_kelvin": temperature_kelvin,
            "n_atoms": ref_trajectory.n_atoms,
        }
        if not dry_run:  # rollout cost provenance (absent for legacy trajectories)
            row |= run_info_row(read_run_info(pred_traj_path))
        rows.append(row | system_metrics)

    df_md = pd.DataFrame(rows).set_index("system")
    if dry_run:
        print(f"Dry run OK for {model_key}: pipeline ran on {system_names}")
        return df_md

    # suffix avoids collisions when parallel single-system jobs share out_dir
    suffix = f"-{'-'.join(system_names)}" if systems is not None else ""
    csv_path = f"{out_dir}/{today}-{model_key}-md-metrics{suffix}.csv.gz"
    df_md.to_csv(csv_path)
    print(f"Per-system MD metrics for {model_key} saved to {csv_path!r}")
    return df_md
