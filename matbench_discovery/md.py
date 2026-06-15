"""Molecular dynamics rollout helpers for the MD benchmark task.

Runs NVT simulations with MLIP calculators using the protocol from the CFPMD-26
benchmark paper: Nose-Hoover chain thermostat, 0.25 fs timestep, 25 fs thermostat
time scale, frames recorded every 10 steps for 20 ps.
"""

import contextlib
import os
import re
import time
import zipfile
from glob import glob

import ase
import ase.io
import numpy as np
import pandas as pd
from ase import Atoms, units
from ase.calculators.calculator import Calculator, PropertyNotImplementedError
from ase.calculators.singlepoint import SinglePointCalculator
from ase.md.nose_hoover_chain import NoseHooverChainNVT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from tqdm import tqdm

from matbench_discovery.metrics import md as md_metrics
from matbench_discovery.trajectory import Trajectory


def extract_temperature(name: str) -> float | None:
    """Temperature in Kelvin parsed from a '<number>K' token in a system name like
    'bulkAu_1500K_Kapil', or None if no such token is found.
    """
    if match := re.search(r"(?<![a-zA-Z0-9.])(\d+(?:\.\d+)?)K(?![a-zA-Z0-9])", name):
        return float(match[1])
    return None


def load_frame_intervals(csv_path: str) -> dict[str, float]:
    """Map '<system>_<temp>K' keys to saved-frame intervals in fs from a CFPMD-style
    settings CSV with System, temperature and dt columns, where dt is the time
    between saved reference frames. Match reference trajectory directories like
    'bulkAg_600K_Kapil' by prefix against the returned keys.
    """
    return {
        f"{row['System']}_{float(row['temperature']):g}K": float(row["dt"])
        for _, row in pd.read_csv(csv_path).iterrows()
    }


def resolve_frame_interval(
    system_name: str, frame_intervals: dict[str, float]
) -> float | None:
    """Return the saved-frame interval for a trajectory directory name.

    CFPMD trajectory directories append suffixes like '_Kapil' or '-Artrith_VASP' to
    settings keys of the form '<system>_<temp>K' (the delimiter after the key is '_'
    for most systems but '-' for a few). Pick the longest delimiter-aware match so
    systems whose names contain another '<...>K' token don't resolve to a shorter
    prefix.
    """
    matches = {
        key: frame_interval
        for key, frame_interval in frame_intervals.items()
        if system_name == key
        or (system_name.startswith(key) and system_name[len(key)] in "_-")
    }
    return matches[max(matches, key=len)] if matches else None


def default_md_reference_paths() -> tuple[str, str]:
    """Reference trajectory dir and settings CSV of the auto-downloaded (and
    auto-extracted) CFPMD-26 dataset.
    """
    from matbench_discovery.enums import DataFiles

    data_root = DataFiles.aimd_reference_md_trajectories.path
    return (
        f"{data_root}/reference_AIMD_trajectories",
        f"{data_root}/reference_AIMD_timestep_and_stride.csv",
    )


def read_trajectory(file_path: str, *, index: str = ":") -> list[Atoms]:
    """Frames of an ASE-readable (optionally compressed) trajectory file. ``index``
    is an ASE slice string, e.g. ':' for all frames or ':64' for the first 64.
    """
    frames = ase.io.read(file_path, index=index)
    return [frames] if isinstance(frames, Atoms) else list(frames)


def run_nvt_md(
    atoms: Atoms,
    calculator: Calculator,
    *,
    temperature_kelvin: float,
    n_steps: int = 80_000,
    time_step_fs: float = 0.25,
    record_interval: int = 10,
    thermostat_time_scale_fs: float = 25,
    seed: int = 0,
    progress_interval: int = 1_000,
) -> list[Atoms]:
    """Run one NVT MD simulation and return recorded frames.

    Velocities are initialized from a Maxwell-Boltzmann distribution. Each recorded
    frame carries energy, forces and (when the calculator supports it) stress via an
    attached SinglePointCalculator, plus velocities, so trajectories written to
    extxyz retain everything needed for RDF/VDOS/pressure evaluation.

    Args:
        atoms: Initial structure (not modified).
        calculator: ASE calculator providing energies, forces and ideally stress.
        temperature_kelvin: Target temperature in Kelvin.
        n_steps: Number of MD steps. Default 80,000 = 20 ps at 0.25 fs.
        time_step_fs: Integration time step in femtoseconds.
        record_interval: Record a frame every this many steps.
        thermostat_time_scale_fs: Nose-Hoover damping time scale in femtoseconds.
        seed: Seed for the Maxwell-Boltzmann velocity initialization.
        progress_interval: Print steps/s and ETA every this many steps (useful for
            monitoring long rollouts in batch job logs). 0 disables progress logs.

    Returns:
        list[Atoms]: Recorded frames including the initial one, i.e.
            n_steps // record_interval + 1 frames.
    """
    atoms = atoms.copy()
    atoms.calc = calculator
    MaxwellBoltzmannDistribution(
        atoms, temperature_K=temperature_kelvin, rng=np.random.default_rng(seed)
    )
    dynamics = _init_nvt_dynamics(
        atoms,
        temperature_kelvin=temperature_kelvin,
        time_step_fs=time_step_fs,
        thermostat_time_scale_fs=thermostat_time_scale_fs,
    )

    frames: list[Atoms] = []
    # ASE calls observers at step 0 and every record_interval steps thereafter
    dynamics.attach(
        lambda: frames.append(_snapshot_frame(atoms, dynamics.nsteps)),
        interval=record_interval,
    )
    _attach_progress_logger(dynamics, n_steps=n_steps, interval=progress_interval)

    dynamics.run(n_steps)
    return frames


def _init_nvt_dynamics(
    atoms: Atoms,
    *,
    temperature_kelvin: float,
    time_step_fs: float,
    thermostat_time_scale_fs: float,
) -> NoseHooverChainNVT:
    """Construct the Nose-Hoover chain NVT integrator (velocities set by caller)."""
    return NoseHooverChainNVT(
        atoms,
        timestep=time_step_fs * units.fs,
        temperature_K=temperature_kelvin,
        tdamp=thermostat_time_scale_fs * units.fs,
    )


def _snapshot_frame(atoms: Atoms, md_step: int) -> Atoms:
    """Copy of the current state with energy/forces/(stress) + velocities attached."""
    results = {"energy": atoms.get_potential_energy(), "forces": atoms.get_forces()}
    with contextlib.suppress(PropertyNotImplementedError):
        results["stress"] = atoms.get_stress(voigt=True)
    frame = atoms.copy()
    frame.info["md_step"] = md_step
    frame.calc = SinglePointCalculator(frame, **results)
    return frame


def _attach_progress_logger(
    dynamics: NoseHooverChainNVT, *, n_steps: int, interval: int
) -> None:
    """Attach a throughput + ETA logger every ``interval`` steps (0 disables)."""
    if interval <= 0:
        return
    start_time = time.perf_counter()

    def log_progress() -> None:
        if (steps_done := dynamics.nsteps) == 0:
            return
        steps_per_sec = steps_done / (time.perf_counter() - start_time)
        eta_min = (n_steps - steps_done) / steps_per_sec / 60
        print(
            f"MD step {steps_done:,}/{n_steps:,}, {steps_per_sec:.1f} steps/s, "
            f"ETA {eta_min:.1f} min",
            flush=True,
        )

    dynamics.attach(log_progress, interval=interval)


# === crash-safe checkpointing for resumable rollouts ===
# NoseHooverChainNVT is deterministic after the seeded velocity init (no per-step RNG),
# so its full state is these five fields; restoring them reproduces an uninterrupted
# run to ~1e-14. Schema-tagged + ASE-version-gated so a renamed attr or upgrade fails
# closed (recompute) rather than silently corrupting trajectories.
CHECKPOINT_SCHEMA = 1


def _nvt_state_attrs(
    dynamics: NoseHooverChainNVT,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The integrator's mutable state arrays; raises if ASE renamed the attributes."""
    try:
        thermostat = dynamics._thermostat  # noqa: SLF001
        state = dynamics._q, dynamics._p, thermostat._eta, thermostat._p_eta  # noqa: SLF001
    except AttributeError as exc:
        raise RuntimeError(
            f"ASE {ase.__version__} NoseHooverChainNVT state attributes changed "
            f"({exc}); MD checkpoint/resume needs updating"
        ) from exc
    return state


def _save_checkpoint(
    ckpt_path: str,
    dynamics: NoseHooverChainNVT,
    *,
    n_steps: int,
    n_frames: int,
    record_interval: int,
    time_step_fs: float,
    temperature_kelvin: float,
    thermostat_time_scale_fs: float,
    seed: int,
) -> None:
    """Atomically write the integrator state + run metadata to ``ckpt_path``."""
    pos, mom, eta, p_eta = _nvt_state_attrs(dynamics)
    tmp_path = f"{ckpt_path}.tmp"
    with open(tmp_path, mode="wb") as file:
        np.savez(
            file,
            schema=CHECKPOINT_SCHEMA,
            ase_version=ase.__version__,
            nsteps=dynamics.nsteps,
            n_steps=n_steps,
            n_frames=n_frames,
            q=pos,
            p=mom,
            eta=eta,
            p_eta=p_eta,
            record_interval=record_interval,
            time_step_fs=time_step_fs,
            temperature_kelvin=temperature_kelvin,
            thermostat_time_scale_fs=thermostat_time_scale_fs,
            seed=seed,
            symbols=np.array(dynamics.atoms.get_chemical_symbols()),
            cell=dynamics.atoms.cell.array,
            pbc=dynamics.atoms.pbc,
        )
    os.replace(tmp_path, ckpt_path)


def _load_checkpoint(
    ckpt_path: str,
    *,
    atoms: Atoms,
    n_steps: int,
    record_interval: int,
    time_step_fs: float,
    temperature_kelvin: float,
    thermostat_time_scale_fs: float,
    seed: int,
) -> dict[str, np.ndarray] | None:
    """Load a checkpoint if present, readable, schema/ASE/parameter-matched; else None
    (so the caller recomputes from scratch rather than resuming a stale/foreign state).
    """
    if not os.path.isfile(ckpt_path):
        return None
    try:
        with np.load(ckpt_path, allow_pickle=False) as data:
            ckpt = {key: data[key] for key in data.files}
    except (OSError, ValueError, EOFError, zipfile.BadZipFile) as exc:
        print(f"  ignoring unreadable checkpoint {ckpt_path}: {exc}")
        return None

    target_frames = n_steps // record_interval + 1
    try:
        nsteps = int(ckpt["nsteps"])
        n_frames = int(ckpt["n_frames"])
        consistent = (
            int(ckpt["schema"]) == CHECKPOINT_SCHEMA
            and str(ckpt["ase_version"]) == ase.__version__
            and int(ckpt["record_interval"]) == record_interval
            and float(ckpt["time_step_fs"]) == time_step_fs
            and float(ckpt["temperature_kelvin"]) == temperature_kelvin
            and float(ckpt["thermostat_time_scale_fs"]) == thermostat_time_scale_fs
            and int(ckpt["seed"]) == seed
            and nsteps <= n_steps
            and n_frames <= target_frames
            and nsteps == (n_frames - 1) * record_interval
            and list(ckpt["symbols"]) == atoms.get_chemical_symbols()
            and np.allclose(ckpt["cell"], atoms.cell.array)
            and bool((ckpt["pbc"] == atoms.pbc).all())
        )
    except (KeyError, TypeError, ValueError):
        consistent = False
    if not consistent:
        print(f"  ignoring checkpoint {ckpt_path}: schema/version/parameters mismatch")
        return None
    return ckpt


def _restore_nvt_state(
    dynamics: NoseHooverChainNVT, ckpt: dict[str, np.ndarray]
) -> None:
    """Restore integrator + thermostat + atoms state from a loaded checkpoint."""
    _nvt_state_attrs(dynamics)  # validate the attributes exist before assigning
    dynamics._q = ckpt["q"]  # noqa: SLF001
    dynamics._p = ckpt["p"]  # noqa: SLF001
    dynamics._thermostat._eta = ckpt["eta"]  # noqa: SLF001
    dynamics._thermostat._p_eta = ckpt["p_eta"]  # noqa: SLF001
    dynamics.nsteps = int(ckpt["nsteps"])
    dynamics.atoms.set_positions(ckpt["q"])
    dynamics.atoms.set_momenta(ckpt["p"])


def run_nvt_md_resumable(
    atoms: Atoms,
    calculator: Calculator,
    *,
    out_path: str,
    temperature_kelvin: float,
    n_steps: int = 80_000,
    time_step_fs: float = 0.25,
    record_interval: int = 10,
    thermostat_time_scale_fs: float = 25,
    seed: int = 0,
    progress_interval: int = 1_000,
    checkpoint_every_n_frames: int = 200,
) -> list[Atoms]:
    """Run an NVT rollout that survives interruption (timeouts/preemption).

    Frames stream to ``{out_path}.part`` and the integrator state is checkpointed to
    ``{out_path}.ckpt.npz`` every ``checkpoint_every_n_frames`` frames. A rerun resumes
    from the last checkpoint (truncating any un-checkpointed tail), continuing the same
    deterministic trajectory. Only on completion is ``.part`` atomically promoted to
    ``out_path`` and the checkpoint removed. Returns the full trajectory.
    """
    # .part.extxyz (not .part) so ASE infers the extxyz format on read/write/append
    part_path = f"{out_path.removesuffix('.extxyz')}.part.extxyz"
    ckpt_path = f"{out_path}.ckpt.npz"
    target_frames = n_steps // record_interval + 1
    if checkpoint_every_n_frames < 1:
        raise ValueError(
            f"{checkpoint_every_n_frames=} must be positive for MD checkpointing"
        )

    atoms = atoms.copy()
    atoms.calc = calculator
    MaxwellBoltzmannDistribution(
        atoms, temperature_K=temperature_kelvin, rng=np.random.default_rng(seed)
    )
    dynamics = _init_nvt_dynamics(
        atoms,
        temperature_kelvin=temperature_kelvin,
        time_step_fs=time_step_fs,
        thermostat_time_scale_fs=thermostat_time_scale_fs,
    )

    ckpt = _load_checkpoint(
        ckpt_path,
        atoms=atoms,
        n_steps=n_steps,
        record_interval=record_interval,
        time_step_fs=time_step_fs,
        temperature_kelvin=temperature_kelvin,
        thermostat_time_scale_fs=thermostat_time_scale_fs,
        seed=seed,
    )
    n_done = 0
    if ckpt is not None and os.path.isfile(part_path):
        n_ckpt = int(ckpt["n_frames"])
        # keep only checkpointed frames; a kill may have appended an un-checkpointed
        # (possibly malformed) tail, so read just the first n_ckpt and rewrite cleanly
        try:
            kept = read_trajectory(part_path, index=f":{n_ckpt}")
        except (OSError, ValueError, IndexError) as exc:
            print(f"  resume failed reading {part_path} ({exc}); recomputing")
            kept = []
        if len(kept) == n_ckpt:
            _restore_nvt_state(dynamics, ckpt)
            ase.io.write(part_path, kept)  # truncate to the checkpointed frames
            n_done = n_ckpt
            print(f"  resuming {out_path} from frame {n_done}/{target_frames}")

    if n_done == 0:  # fresh start: clear any stale partial files
        for stale in (part_path, ckpt_path):
            with contextlib.suppress(FileNotFoundError):
                os.remove(stale)

    frames_written = n_done

    def record_and_checkpoint() -> None:
        nonlocal frames_written
        ase.io.write(part_path, _snapshot_frame(atoms, dynamics.nsteps), append=True)
        frames_written += 1
        if frames_written % checkpoint_every_n_frames == 0:
            _save_checkpoint(
                ckpt_path,
                dynamics,
                n_steps=n_steps,
                n_frames=frames_written,
                record_interval=record_interval,
                time_step_fs=time_step_fs,
                temperature_kelvin=temperature_kelvin,
                thermostat_time_scale_fs=thermostat_time_scale_fs,
                seed=seed,
            )

    dynamics.attach(record_and_checkpoint, interval=record_interval)
    _attach_progress_logger(dynamics, n_steps=n_steps, interval=progress_interval)

    dynamics.run(n_steps - dynamics.nsteps)  # remaining steps (0 nsteps if fresh)

    frames = read_trajectory(part_path)
    if len(frames) != target_frames:
        raise RuntimeError(
            f"{out_path} rollout produced {len(frames)} frames != {target_frames}"
        )
    os.replace(part_path, out_path)  # atomic promote only once complete
    with contextlib.suppress(FileNotFoundError):
        os.remove(ckpt_path)
    return frames


def find_reference_trajectory(ref_dir: str, system_name: str) -> str:
    """Path to the single reference trajectory in ref_dir/system_name, matching
    traj.extxyz and ASE-transparent .gz/.xz/.bz2 compressed variants.
    """
    traj_files = sorted(glob(f"{ref_dir}/{system_name}/traj.*xyz*"))
    if len(traj_files) != 1:
        raise ValueError(
            f"Expected 1 reference trajectory in {ref_dir}/{system_name}, "
            f"got {traj_files}"
        )
    return traj_files[0]


def read_reference_trajectory(
    ref_dir: str, system_name: str, *, max_frames: int | None = None
) -> Trajectory:
    """Reference trajectory as a Trajectory, reading the fast HDF5 artifact
    (``<system>/traj.h5``) when present and falling back to parsing extxyz otherwise.

    The HDF5 path only reads/decompresses the requested frames, so huge references
    load in milliseconds; the extxyz fallback (ASE's slow parser) keeps the pipeline
    working on un-converted datasets. ``max_frames`` caps the leading frames read
    (used by dry runs to stay fast).
    """
    h5_path = f"{ref_dir}/{system_name}/traj.h5"
    if os.path.isfile(h5_path):
        frames = slice(0, max_frames) if max_frames is not None else slice(None)
        return Trajectory.read_hdf5(h5_path, frames=frames)
    ext_path = find_reference_trajectory(ref_dir, system_name)
    index = f":{max_frames}" if max_frames is not None else ":"
    atoms_frames = read_trajectory(ext_path, index=index)
    if not atoms_frames:
        raise ValueError(
            f"Reference trajectory {ext_path!r} for {system_name!r} empty or corrupted"
        )
    return Trajectory.from_ase(atoms_frames)


def validate_pred_trajectory(
    path: str, *, expected_frames: int, ref_atoms: Atoms, record_interval: int
) -> list[Atoms] | None:
    """Read and validate a cached predicted rollout, returning its frames if sound or
    None if it should be recomputed.

    A truncated or otherwise inconsistent rollout would silently evaluate as a shorter
    trajectory (evaluate_md_system time-matches to the shorter span), so reuse only
    when the file is complete and consistent: readable, exactly ``expected_frames``
    frames with a strict ``md_step`` sequence 0, record_interval, 2*record_interval,
    ...; matching atom count and chemical symbols, cell and pbc of ``ref_atoms``; and
    all-finite positions and velocities.
    """
    try:
        frames = read_trajectory(path)
    except (OSError, ValueError, IndexError, StopIteration) as exc:
        print(f"  recomputing {path}: unreadable ({exc})")
        return None

    if (n_frames := len(frames)) != expected_frames:
        print(f"  recomputing {path}: {n_frames} frames != expected {expected_frames}")
        return None

    ref_symbols = ref_atoms.get_chemical_symbols()
    for idx, frame in enumerate(frames):
        velocities = frame.get_velocities()
        reason = None
        if int(frame.info.get("md_step", -1)) != idx * record_interval:
            reason = f"md_step {frame.info.get('md_step')} != {idx * record_interval}"
        elif frame.get_chemical_symbols() != ref_symbols:
            reason = "chemical symbols differ from reference"
        elif not np.allclose(frame.cell.array, ref_atoms.cell.array) or bool(
            (frame.pbc != ref_atoms.pbc).any()
        ):
            reason = "cell or pbc differ from reference"
        elif not np.isfinite(frame.positions).all() or (
            velocities is not None and not np.isfinite(velocities).all()
        ):
            reason = "non-finite positions or velocities"
        if reason:
            print(f"  recomputing {path}: frame {idx} {reason}")
            return None

    return frames


def run_md_benchmark(
    *,
    calculator: Calculator,
    model_key: str,
    out_dir: str,
    ref_dir: str | None = None,
    settings_csv: str | None = None,
    ref_frame_interval_fs: float | None = None,
    systems: "list[str] | None" = None,
    n_steps: int = 80_000,
    time_step_fs: float = 0.25,
    record_interval: int = 10,
    seed: int = 0,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Run NVT rollouts and compute MD metrics for one model across reference systems.

    For each system directory under ref_dir, this rolls out an NVT trajectory from
    the reference initial structure (reusing an existing rollout if present), then
    compares it to the reference via energy/force RMSE, RDF, VDOS and pressure. Writes
    one gzipped per-system metrics CSV and returns the per-system DataFrame.

    Args:
        calculator: ASE calculator for the model under test.
        model_key: Model enum name/key, used in output filenames.
        out_dir: Directory for rollout trajectories and the metrics CSV.
        ref_dir: Reference trajectory dir. Defaults to the auto-downloaded CFPMD-26
            dataset (and settings_csv defaults alongside it).
        settings_csv: Per-system frame-interval CSV. Mutually informs ref dt.
        ref_frame_interval_fs: Constant reference frame interval (overrides CSV).
        systems: Subset of system dir names to run. Defaults to all under ref_dir.
        n_steps: MD steps per rollout (80,000 = 20 ps at 0.25 fs).
        time_step_fs: MD integration time step.
        record_interval: Record a frame every this many MD steps.
        seed: Maxwell-Boltzmann velocity seed.
        dry_run: Smoke test: one system, a handful of MD steps, and a capped reference
            slice, to verify the model + pipeline run error-free in seconds without
            writing outputs. Returned metrics are not meaningful.

    Returns:
        pd.DataFrame: One row per system indexed by system name.
    """
    if ref_dir is None:
        ref_dir, default_csv = default_md_reference_paths()
        settings_csv = settings_csv or default_csv
    elif settings_csv is None and ref_frame_interval_fs is None:
        raise ValueError(
            "Pass settings_csv or ref_frame_interval_fs for custom ref_dir"
        )

    frame_intervals = load_frame_intervals(settings_csv) if settings_csv else {}
    system_dirs = sorted(
        entry.name
        for entry in os.scandir(ref_dir)
        if entry.is_dir() and (not systems or entry.name in systems)
    )
    if not system_dirs:  # clearer than the downstream empty-DataFrame KeyError
        raise ValueError(
            f"No system directories found under {ref_dir!r}"
            + (f" matching {systems=}" if systems else "")
        )
    if dry_run:  # one system, a few steps, capped reference slice
        system_dirs = system_dirs[:1]
        n_steps = record_interval * 8
    else:
        os.makedirs(out_dir, exist_ok=True)  # for rollout trajectories and metrics CSV

    pred_time_step_fs = time_step_fs * record_interval
    rows: list[dict[str, float | str]] = []
    for system_name in tqdm(system_dirs, desc=f"MD systems ({model_key})"):
        temperature_kelvin = extract_temperature(system_name)
        if temperature_kelvin is None:
            raise ValueError(f"Could not parse temperature from {system_name!r}")
        ref_dt_fs = ref_frame_interval_fs or resolve_frame_interval(
            system_name, frame_intervals
        )
        if ref_dt_fs is None:
            raise ValueError(f"No frame interval for {system_name!r} in settings CSV")

        # dry run reads only the first 64 reference frames so the single-point eval
        # stays fast; the HDF5 reader only decompresses the frames it returns
        ref_trajectory = read_reference_trajectory(
            ref_dir, system_name, max_frames=64 if dry_run else None
        )
        initial_atoms = ref_trajectory.frame_as_atoms(0)

        pred_traj_path = f"{out_dir}/{system_name}-nvt-{model_key}.extxyz"
        target_frames = n_steps // record_interval + 1
        pred_frames = None  # list[Atoms] from reuse, in-memory run, or resumable run
        if not dry_run and os.path.isfile(pred_traj_path):
            # only reuse a complete, consistent rollout; else recompute (a truncated
            # file would silently evaluate as a shorter trajectory)
            pred_frames = validate_pred_trajectory(
                pred_traj_path,
                expected_frames=target_frames,
                ref_atoms=initial_atoms,
                record_interval=record_interval,
            )
        if pred_frames is None and dry_run:  # in-memory, no checkpoint files
            pred_frames = run_nvt_md(
                initial_atoms,
                calculator,
                temperature_kelvin=temperature_kelvin,
                n_steps=n_steps,
                time_step_fs=time_step_fs,
                record_interval=record_interval,
                seed=seed,
            )
        elif pred_frames is None:  # crash-safe resumable rollout to pred_traj_path
            pred_frames = run_nvt_md_resumable(
                initial_atoms,
                calculator,
                out_path=pred_traj_path,
                temperature_kelvin=temperature_kelvin,
                n_steps=n_steps,
                time_step_fs=time_step_fs,
                record_interval=record_interval,
                seed=seed,
            )

        system_metrics = md_metrics.evaluate_md_system(
            ref_trajectory,
            pred_frames,  # evaluate_md_system coerces the Atoms list to a Trajectory
            ref_time_step_fs=ref_dt_fs,
            pred_time_step_fs=pred_time_step_fs,
            calculator=calculator,
        )
        rows.append(
            {"system": system_name, "temperature_kelvin": temperature_kelvin}
            | system_metrics
        )

    df_md = pd.DataFrame(rows).set_index("system")
    if dry_run:
        print(f"Dry run OK for {model_key}: pipeline ran on {system_dirs}")
        return df_md

    # suffix avoids collisions when parallel single-system jobs share out_dir
    suffix = f"-{'-'.join(systems)}" if systems else ""
    csv_path = f"{out_dir}/{model_key}-md-metrics{suffix}.csv.gz"
    df_md.to_csv(csv_path)
    print(f"Per-system MD metrics for {model_key} saved to {csv_path!r}")
    return df_md
