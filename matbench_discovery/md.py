"""Molecular dynamics rollout helpers for the MD benchmark task.

Runs NVT simulations with MLIP calculators using the protocol from the CFPMD-26
benchmark paper: Nose-Hoover chain thermostat, 0.25 fs timestep, 25 fs thermostat
time scale, frames recorded every 10 steps for 20 ps.
"""

import contextlib
import os
import zipfile
from collections.abc import Iterable
from dataclasses import dataclass

import ase
import ase.io
import h5py
import numpy as np
import pandas as pd
from ase import Atoms, units
from ase.calculators.calculator import Calculator, PropertyNotImplementedError
from ase.calculators.singlepoint import SinglePointCalculator
from ase.md.nose_hoover_chain import NoseHooverChainNVT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.metrics import md as md_metrics
from matbench_discovery.trajectory import Trajectory

# === single-file multi-system reference dataset ===
# CFPMD-26 references ship as one HDF5 with a group per system, each carrying the
# trajectory arrays (guarded by Trajectory's per-group TRAJECTORY_SCHEMA) plus its
# saved-frame interval (dt_fs) and target temperature as group attrs. Self-describing
# -> no separate settings CSV or system-name matching.


def default_md_reference_path() -> str:
    """Path to the auto-downloaded CFPMD-26 reference HDF5 (one group per system)."""
    from matbench_discovery.enums import DataFiles

    return DataFiles.aimd_reference_md_trajectories.path


def write_reference_h5(
    path: str, entries: Iterable[tuple[str, Trajectory, float, float]]
) -> None:
    """Atomically write a multi-system reference HDF5: one group per system holding the
    trajectory arrays plus ``dt_fs`` (saved-frame interval) and ``temperature_kelvin``
    attrs. ``entries`` yields (system_name, trajectory, dt_fs, temperature_kelvin); it
    may be a lazy generator so the converter parses one trajectory at a time and writes
    it straight to disk, bounding peak memory to the single largest system.
    """
    if dir_name := os.path.dirname(path):
        os.makedirs(dir_name, exist_ok=True)
    tmp_path = f"{path}.tmp"
    n_written = 0
    with h5py.File(tmp_path, "w") as file:
        for system_name, trajectory, dt_fs, temperature_kelvin in entries:
            group = file.create_group(system_name)
            trajectory.write_to_h5_group(group)
            group.attrs["dt_fs"] = float(dt_fs)
            group.attrs["temperature_kelvin"] = float(temperature_kelvin)
            n_written += 1
    if n_written == 0:  # nothing yielded: drop the empty tmp file and fail loud
        os.remove(tmp_path)
        raise ValueError("Cannot write a reference file with no systems")
    os.replace(tmp_path, path)


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
    progress_label: str = "MD rollout",
) -> list[Atoms]:
    """Run one NVT MD simulation and return recorded frames.

    Velocities are initialized from a Maxwell-Boltzmann distribution. Each recorded
    frame carries energy, forces and (when the calculator supports it) stress via an
    attached SinglePointCalculator, plus velocities, so trajectories written to
    extxyz retain everything needed for RDF/vDOS/pressure evaluation.

    Args:
        atoms: Initial structure (not modified).
        calculator: ASE calculator providing energies, forces and ideally stress.
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
        list[Atoms]: Recorded frames including the initial one, i.e.
            n_steps // record_interval + 1 frames.
    """
    atoms, dynamics = _thermalized_nvt_dynamics(
        atoms,
        calculator,
        temperature_kelvin=temperature_kelvin,
        time_step_fs=time_step_fs,
        thermostat_time_scale_fs=thermostat_time_scale_fs,
        seed=seed,
    )

    frames: list[Atoms] = []
    # ASE calls observers at step 0 and every record_interval steps thereafter
    dynamics.attach(
        lambda: frames.append(_snapshot_frame(atoms, dynamics.nsteps)),
        interval=record_interval,
    )
    _run_dynamics_with_progress(
        dynamics,
        n_steps=n_steps,
        interval=progress_interval,
        label=progress_label,
    )

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


def _thermalized_nvt_dynamics(
    atoms: Atoms,
    calculator: Calculator,
    *,
    temperature_kelvin: float,
    time_step_fs: float,
    thermostat_time_scale_fs: float,
    seed: int,
) -> tuple[Atoms, NoseHooverChainNVT]:
    """Copy ``atoms``, attach the calculator, seed Maxwell-Boltzmann velocities and
    build the NVT integrator. Single deterministic setup shared by the in-memory and
    resumable rollouts so both thermalize identically for a given seed.
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
    return atoms, dynamics


def _snapshot_frame(atoms: Atoms, md_step: int) -> Atoms:
    """Copy of the current state with energy/forces/(stress) + velocities attached."""
    results = {"energy": atoms.get_potential_energy(), "forces": atoms.get_forces()}
    with contextlib.suppress(PropertyNotImplementedError):
        results["stress"] = atoms.get_stress(voigt=True)
    frame = atoms.copy()
    frame.info["md_step"] = md_step
    frame.calc = SinglePointCalculator(frame, **results)
    return frame


def _run_dynamics_with_progress(
    dynamics: NoseHooverChainNVT, *, n_steps: int, interval: int, label: str
) -> None:
    """Run remaining MD steps with a tqdm progress bar (0 interval disables it)."""
    remaining_steps = n_steps - dynamics.nsteps
    if interval <= 0:
        dynamics.run(remaining_steps)
        return

    with tqdm(
        total=n_steps,
        initial=dynamics.nsteps,
        desc=label,
        unit="step",
        mininterval=5,
        leave=True,
    ) as pbar:

        def update_progress() -> None:
            if (step_delta := dynamics.nsteps - pbar.n) > 0:
                pbar.update(step_delta)

        dynamics.attach(update_progress, interval=interval)
        dynamics.run(remaining_steps)
        update_progress()


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
    params: NvtParams,
    n_steps: int,
    n_frames: int,
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
            record_interval=params.record_interval,
            time_step_fs=params.time_step_fs,
            temperature_kelvin=params.temperature_kelvin,
            thermostat_time_scale_fs=params.thermostat_time_scale_fs,
            seed=params.seed,
            symbols=np.array(dynamics.atoms.get_chemical_symbols()),
            cell=dynamics.atoms.cell.array,
            pbc=dynamics.atoms.pbc,
        )
    os.replace(tmp_path, ckpt_path)


def _load_checkpoint(
    ckpt_path: str, *, atoms: Atoms, params: NvtParams, n_steps: int
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

    target_frames = n_steps // params.record_interval + 1
    try:
        nsteps = int(ckpt["nsteps"])
        n_frames = int(ckpt["n_frames"])
        stored = NvtParams(
            temperature_kelvin=float(ckpt["temperature_kelvin"]),
            time_step_fs=float(ckpt["time_step_fs"]),
            record_interval=int(ckpt["record_interval"]),
            thermostat_time_scale_fs=float(ckpt["thermostat_time_scale_fs"]),
            seed=int(ckpt["seed"]),
        )
        consistent = (
            int(ckpt["schema"]) == CHECKPOINT_SCHEMA
            and str(ckpt["ase_version"]) == ase.__version__
            and stored == params
            and nsteps <= n_steps
            and n_frames <= target_frames
            and nsteps == (n_frames - 1) * params.record_interval
            and list(ckpt["symbols"]) == atoms.get_chemical_symbols()
            and np.allclose(ckpt["cell"], atoms.cell.array)
            and np.array_equal(ckpt["pbc"], atoms.pbc)
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
    progress_label: str = "MD rollout",
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

    atoms, dynamics = _thermalized_nvt_dynamics(
        atoms,
        calculator,
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

    ckpt = _load_checkpoint(ckpt_path, atoms=atoms, params=params, n_steps=n_steps)
    frames_written = 0
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
            frames_written = n_ckpt
            print(f"  resuming {out_path} from frame {frames_written}/{target_frames}")

    if frames_written == 0:  # fresh start: clear any stale partial files
        for stale in (part_path, ckpt_path):
            with contextlib.suppress(FileNotFoundError):
                os.remove(stale)

    def record_and_checkpoint() -> None:
        nonlocal frames_written
        ase.io.write(part_path, _snapshot_frame(atoms, dynamics.nsteps), append=True)
        frames_written += 1
        if frames_written % checkpoint_every_n_frames == 0:
            _save_checkpoint(
                ckpt_path,
                dynamics,
                params=params,
                n_steps=n_steps,
                n_frames=frames_written,
            )

    dynamics.attach(record_and_checkpoint, interval=record_interval)
    _run_dynamics_with_progress(
        dynamics,
        n_steps=n_steps,
        interval=progress_interval,
        label=progress_label,
    )

    frames = read_trajectory(part_path)
    if len(frames) != target_frames:
        raise RuntimeError(
            f"{out_path} rollout produced {len(frames)} frames != {target_frames}"
        )
    os.replace(part_path, out_path)  # atomic promote only once complete
    with contextlib.suppress(FileNotFoundError):
        os.remove(ckpt_path)
    return frames


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
        elif not np.allclose(frame.cell.array, ref_atoms.cell.array) or np.any(
            frame.pbc != ref_atoms.pbc
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
    ref_file: str | None = None,
    systems: list[str] | None = None,
    n_steps: int = 80_000,
    time_step_fs: float = 0.25,
    record_interval: int = 10,
    seed: int = 0,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Run NVT rollouts and compute MD metrics for one model across reference systems.

    For each system in the reference HDF5, this rolls out an NVT trajectory from the
    reference initial structure (reusing an existing rollout if present), then compares
    it to the reference via energy/force RMSE, RDF, ADF, vDOS and pressure. Writes one
    gzipped per-system metrics CSV and returns the per-system DataFrame.

    Args:
        calculator: ASE calculator for the model under test.
        model_key: Model enum name/key, used in output filenames.
        out_dir: Directory for rollout trajectories and the metrics CSV.
        ref_file: Reference HDF5 (one group per system, each carrying dt_fs and
            temperature_kelvin attrs). Defaults to the auto-downloaded CFPMD-26 set.
        systems: Subset of system names to run. Defaults to all in ref_file.
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
    if ref_file is None:
        ref_file = default_md_reference_path()

    all_systems = list_reference_systems(ref_file)
    system_names = (
        all_systems
        if systems is None
        else [name for name in all_systems if name in systems]
    )
    if not system_names:  # clearer than the downstream empty-DataFrame KeyError
        raise ValueError(
            f"No systems found in reference {ref_file!r}"
            + (f" matching {systems=}" if systems else "")
        )
    if dry_run:  # one system, a few steps, capped reference slice
        system_names = system_names[:1]
        # enough recorded frames that time-matching keeps >=4 even for the coarsest
        # reference cadence (dt_fs up to 2 fs vs the 2.5 fs prediction cadence)
        n_steps = record_interval * 20
    else:
        os.makedirs(out_dir, exist_ok=True)  # for rollout trajectories and metrics CSV

    pred_time_step_fs = time_step_fs * record_interval
    rows: list[dict[str, float | str]] = []
    for system_name in tqdm(system_names, desc=f"MD systems ({model_key})"):
        progress_label = f"{model_key} {system_name}"
        # dry run reads only the first 64 reference frames so the single-point eval
        # stays fast; the HDF5 reader only decompresses the frames it returns
        ref_trajectory, ref_time_step_fs, temperature_kelvin = (
            read_reference_trajectory(
                ref_file, system_name, max_frames=64 if dry_run else None
            )
        )
        initial_atoms = ref_trajectory.frame_as_atoms(0)

        pred_traj_path = f"{out_dir}/{system_name}-nvt-{model_key}.extxyz"
        target_frames = n_steps // record_interval + 1
        if dry_run:  # in-memory run, no checkpoint/output files
            pred_frames = run_nvt_md(
                initial_atoms,
                calculator,
                temperature_kelvin=temperature_kelvin,
                n_steps=n_steps,
                time_step_fs=time_step_fs,
                record_interval=record_interval,
                seed=seed,
                progress_label=progress_label,
            )
        else:
            # reuse only a complete, consistent rollout; else recompute (a truncated
            # file would silently evaluate as a shorter trajectory)
            pred_frames = None
            if os.path.isfile(pred_traj_path):
                pred_frames = validate_pred_trajectory(
                    pred_traj_path,
                    expected_frames=target_frames,
                    ref_atoms=initial_atoms,
                    record_interval=record_interval,
                )
            if pred_frames is None:  # crash-safe resumable rollout to pred_traj_path
                pred_frames = run_nvt_md_resumable(
                    initial_atoms,
                    calculator,
                    out_path=pred_traj_path,
                    temperature_kelvin=temperature_kelvin,
                    n_steps=n_steps,
                    time_step_fs=time_step_fs,
                    record_interval=record_interval,
                    seed=seed,
                    progress_label=progress_label,
                )

        refeval_path = f"{out_dir}/{system_name}-nvt-{model_key}-refeval.npz"
        ref_predictions: dict[str, np.ndarray] | None = None
        if os.path.isfile(refeval_path):
            with np.load(refeval_path, allow_pickle=False) as cached:
                ref_predictions = {key: cached[key] for key in ("e_pred", "force_se")}
            try:
                md_metrics.energy_force_rmse_from_preds(ref_trajectory, ref_predictions)
            except ValueError as exc:
                print(f"  recomputing {refeval_path}: {exc}", flush=True)
                ref_predictions = None
            else:
                print(f"  reusing {refeval_path}", flush=True)

        if ref_predictions is None:
            ref_predictions = md_metrics.predict_on_reference(
                ref_trajectory, calculator
            )
            if not dry_run:
                # e_ref lives in HDF5; the sidecar stores only model predictions.
                np.savez(
                    refeval_path,
                    e_pred=ref_predictions["e_pred"],
                    force_se=ref_predictions["force_se"],
                    n_atoms=ref_trajectory.n_atoms,
                )
        system_metrics = md_metrics.evaluate_md_system(
            ref_trajectory,
            pred_frames,  # evaluate_md_system coerces the Atoms list to a Trajectory
            ref_time_step_fs=ref_time_step_fs,
            pred_time_step_fs=pred_time_step_fs,
            ref_predictions=ref_predictions,
            progress_label=progress_label,
        )
        rows.append(
            {"system": system_name, "temperature_kelvin": temperature_kelvin}
            | system_metrics
        )

    df_md = pd.DataFrame(rows).set_index("system")
    if dry_run:
        print(f"Dry run OK for {model_key}: pipeline ran on {system_names}")
        return df_md

    # suffix avoids collisions when parallel single-system jobs share out_dir
    suffix = f"-{'-'.join(systems)}" if systems else ""
    csv_path = f"{out_dir}/{today}-{model_key}-md-metrics{suffix}.csv.gz"
    df_md.to_csv(csv_path)
    print(f"Per-system MD metrics for {model_key} saved to {csv_path!r}")
    return df_md
