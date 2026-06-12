"""Molecular dynamics rollout helpers for the MD benchmark task.

Runs NVT simulations with MLIP calculators using the protocol from the CFPMD-26
benchmark paper: Nose-Hoover chain thermostat, 0.25 fs timestep, 25 fs thermostat
time scale, frames recorded every 10 steps for 20 ps.
"""

import contextlib
import re
import time

import ase.io
import numpy as np
import pandas as pd
from ase import Atoms, units
from ase.calculators.calculator import Calculator, PropertyNotImplementedError
from ase.calculators.singlepoint import SinglePointCalculator
from ase.md.nose_hoover_chain import NoseHooverChainNVT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


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

    CFPMD trajectory directories append suffixes like '_Kapil' to settings keys of
    the form '<system>_<temp>K'. Pick the longest delimiter-aware match so systems
    whose names contain another '<...>K' token don't resolve to a shorter prefix.
    """
    matches = {
        key: frame_interval
        for key, frame_interval in frame_intervals.items()
        if system_name == key or system_name.startswith(f"{key}_")
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


def read_trajectory(file_path: str) -> list[Atoms]:
    """All frames of an ASE-readable (optionally compressed) trajectory file."""
    frames = ase.io.read(file_path, index=":")
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
