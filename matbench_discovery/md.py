"""Molecular dynamics rollout helpers for the MD benchmark task.

Runs NVT simulations with MLIP calculators using the protocol from the CFPMD-26
benchmark paper: Nose-Hoover chain thermostat, 0.25 fs timestep, 25 fs thermostat
time scale, frames recorded every 10 steps for 20 ps.
"""

import contextlib
import os
import re
import time
from glob import glob

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
        f"{row['System']}_{int(row['temperature'])}K": float(row["dt"])
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
    dynamics = NoseHooverChainNVT(
        atoms,
        timestep=time_step_fs * units.fs,
        temperature_K=temperature_kelvin,
        tdamp=thermostat_time_scale_fs * units.fs,
    )

    frames: list[Atoms] = []

    def record_frame() -> None:
        """Snapshot the current state with results attached."""
        results = {
            "energy": atoms.get_potential_energy(),
            "forces": atoms.get_forces(),
        }
        with contextlib.suppress(PropertyNotImplementedError):
            results["stress"] = atoms.get_stress(voigt=True)
        frame = atoms.copy()
        frame.info["md_step"] = dynamics.nsteps
        frame.calc = SinglePointCalculator(frame, **results)
        frames.append(frame)

    # ASE calls observers at step 0 and every record_interval steps thereafter
    dynamics.attach(record_frame, interval=record_interval)

    if progress_interval > 0:
        start_time = time.perf_counter()

        def log_progress() -> None:
            """Print MD throughput and remaining-time estimate."""
            if (steps_done := dynamics.nsteps) == 0:
                return
            steps_per_sec = steps_done / (time.perf_counter() - start_time)
            eta_min = (n_steps - steps_done) / steps_per_sec / 60
            print(
                f"MD step {steps_done:,}/{n_steps:,}, {steps_per_sec:.1f} steps/s, "
                f"ETA {eta_min:.1f} min",
                flush=True,
            )

        dynamics.attach(log_progress, interval=progress_interval)

    dynamics.run(n_steps)
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
        # (and the slow xz decompression of multi-GB references) stays fast
        ref_path = find_reference_trajectory(ref_dir, system_name)
        ref_trajectory = read_trajectory(ref_path, index=":64" if dry_run else ":")

        pred_traj_path = f"{out_dir}/{system_name}-nvt-{model_key}.extxyz"
        if not dry_run and os.path.isfile(pred_traj_path):
            pred_trajectory = read_trajectory(pred_traj_path)
        else:
            pred_trajectory = run_nvt_md(
                ref_trajectory[0],
                calculator,
                temperature_kelvin=temperature_kelvin,
                n_steps=n_steps,
                time_step_fs=time_step_fs,
                record_interval=record_interval,
                seed=seed,
            )
            if not dry_run:
                ase.io.write(pred_traj_path, pred_trajectory)

        system_metrics = md_metrics.evaluate_md_system(
            ref_trajectory,
            pred_trajectory,
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
