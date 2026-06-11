"""Metrics for evaluating models on reference molecular dynamics trajectories."""

from __future__ import annotations

import csv
import json
import math
import re
import traceback
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import ase.io
import numpy as np
import pandas as pd
from ruamel.yaml.comments import CommentedMap
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.data import update_yaml_file
from matbench_discovery.md import (
    MdRunConfig,
    find_trajectory_files,
    infer_temperature,
    read_trajectory,
    should_skip_trajectory,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ase import Atoms
    from ase.calculators.calculator import Calculator

    from matbench_discovery.enums import Model

ENERGY_INFO_KEYS = (
    "energy",
    "Energy",
    "dft_energy",
    "DFT_energy",
    "ref_energy",
    "reference_energy",
)
FORCE_ARRAY_KEYS = (
    "forces",
    "force",
    "dft_forces",
    "DFT_forces",
    "ref_forces",
    "reference_forces",
)

RdfMode = Literal["same", "different"]
VdosMode = Literal["same", "different"]
CM_INV_TO_EV = 1.2398419843320026e-4
SPEED_OF_LIGHT_CM_S = 2.9979245899e10
EV_PER_A3_TO_GPA = 160.21766208


@dataclass(frozen=True)
class MdRdfConfig:
    """Configuration for RDF trajectory error metrics."""

    mode: RdfMode = "same"
    nbins: int = 500
    ref_frame_dt_fs: float | None = None
    mlip_frame_dt_fs: float | None = None
    ref_settings_path: str | Path | None = None
    mlip_settings_path: str | Path | None = None
    write_rdf_curves: bool = True


@dataclass(frozen=True)
class MdRdfEvalResult:
    """Per-system RDF error between reference and MLIP MD trajectories."""

    model_key: str
    system_name: str
    ref_file: str
    mlip_file: str | None
    rdf_error: float | None
    n_ref_frames: int
    n_mlip_frames: int
    ref_frame_dt_fs: float | None
    mlip_frame_dt_fs: float | None
    matched_time_fs: float | None
    status: str = "completed"
    error: str | None = None


@dataclass(frozen=True)
class MdVdosConfig:
    """Configuration for source-compatible VDOS error metrics.

    Frame spacings are the native trajectory intervals before the optional
    VDOS stride is applied. The effective VDOS spacing is ``frame_dt * stride``.
    """

    mode: VdosMode = "same"
    ref_frame_dt_fs: float | None = None
    mlip_frame_dt_fs: float | None = None
    ref_stride: int | None = None
    mlip_stride: int | None = None
    ref_padding: int | None = None
    mlip_padding: int | None = None
    ref_settings_path: str | Path | None = None
    mlip_settings_path: str | Path | None = None
    e_min: float | None = None
    e_max: float | None = None
    normalize_area: bool = True
    write_vdos_curves: bool = True


@dataclass(frozen=True)
class MdVdosEvalResult:
    """Per-system normalized VDOS error between MD trajectories."""

    model_key: str
    system_name: str
    ref_file: str
    mlip_file: str | None
    vdos_error: float | None
    n_ref_frames: int
    n_mlip_frames: int
    ref_frame_dt_fs: float | None
    mlip_frame_dt_fs: float | None
    ref_stride: int | None
    mlip_stride: int | None
    ref_padding: int | None
    mlip_padding: int | None
    matched_time_fs: float | None
    status: str = "completed"
    error: str | None = None


@dataclass(frozen=True)
class MdPressureConfig:
    """Configuration for same-simulation-length pressure histogram metrics."""

    bins: int = 80
    ref_frame_dt_fs: float | None = None
    mlip_frame_dt_fs: float | None = None
    ref_settings_path: str | Path | None = None
    mlip_settings_path: str | Path | None = None
    two_d_systems: tuple[str, ...] = ("TiSe2",)
    plane_2d: Literal["xy", "xz", "yz"] = "xy"
    skip_dir_prefixes: tuple[str, ...] = ("Pt111w24H2O_",)
    write_pressure_values: bool = True


@dataclass(frozen=True)
class MdPressureEvalResult:
    """Per-system pressure distribution error between reference and MLIP MD."""

    model_key: str
    system_name: str
    ref_file: str
    mlip_file: str | None
    pressure_error_percent: float | None
    pressure_histogram_l1_area: float | None
    pressure_histogram_distance: float | None
    pressure_ref_histogram_area: float | None
    pressure_mlip_histogram_area: float | None
    pressure_histogram_min_GPa: float | None  # noqa: N815
    pressure_histogram_max_GPa: float | None  # noqa: N815
    bins: int
    n_ref_frames: int
    n_mlip_frames: int
    ref_frame_dt_fs: float | None
    mlip_frame_dt_fs: float | None
    matched_time_fs: float | None
    pressure_mode: str | None = None
    plane_used: str | None = None
    status: str = "completed"
    error: str | None = None


@dataclass(frozen=True)
class MdReferenceEvalResult:
    """Per-frame prediction errors on a reference MD trajectory."""

    model_key: str
    input_file: str
    system_name: str
    frame_idx: int
    n_atoms: int
    temperature_k: float | None
    ref_energy: float | None
    pred_energy: float | None
    energy_error_per_atom: float | None
    force_rmse: float | None
    force_mae: float | None
    force_max_error: float | None
    force_sse: float | None
    n_force_components: int
    status: str = "completed"
    error: str | None = None


def get_reference_energy(atoms: Atoms) -> float | None:
    """Get reference potential energy from ASE atoms metadata or calculator results."""
    for key in ENERGY_INFO_KEYS:
        if key not in atoms.info:
            continue
        value = np.asarray(atoms.info[key], dtype=float)
        if value.size == 1:
            return float(value.reshape(-1)[0])

    try:
        return float(atoms.get_potential_energy())
    except Exception:  # noqa: BLE001 - ASE calculators expose diverse failures
        return None


def get_reference_forces(atoms: Atoms) -> np.ndarray | None:
    """Get reference forces from ASE atoms arrays or calculator results."""
    for key in FORCE_ARRAY_KEYS:
        if key not in atoms.arrays:
            continue
        forces = np.asarray(atoms.arrays[key], dtype=float)
        if forces.shape == (len(atoms), 3):
            return forces

    try:
        forces = np.asarray(atoms.get_forces(), dtype=float)
    except Exception:  # noqa: BLE001 - ASE calculators expose diverse failures
        return None

    if forces.shape != (len(atoms), 3):
        return None
    return forces


def evaluate_reference_frames(
    frames: Sequence[Atoms],
    *,
    calculator: Calculator,
    model_key: str,
    input_file: str | Path,
    system_name: str,
    progress_bar: tqdm | None = None,
) -> pd.DataFrame:
    """Evaluate model energy and forces on reference MD frames.

    Energies are compared per atom in eV/atom. Force errors are computed over all
    Cartesian force components in eV/Angstrom.
    """
    rows: list[dict[str, Any]] = []

    for frame_idx, atoms in enumerate(frames):
        if progress_bar is not None:
            progress_bar.set_postfix_str(
                f"{system_name} frame {frame_idx + 1}/{len(frames)}",
                refresh=False,
            )

        ref_energy = get_reference_energy(atoms)
        ref_forces = get_reference_forces(atoms)
        temperature_k = infer_temperature(input_file, atoms=atoms)

        if ref_energy is None and ref_forces is None:
            rows.append(
                {
                    "model_key": model_key,
                    "input_file": str(input_file),
                    "system_name": system_name,
                    "frame_idx": frame_idx,
                    "n_atoms": len(atoms),
                    "temperature_k": temperature_k,
                    "status": "missing_reference",
                    "error": "No reference energy or forces found",
                }
            )
            if progress_bar is not None:
                progress_bar.update()
            continue

        pred_atoms = atoms.copy()
        pred_atoms.calc = calculator

        pred_energy: float | None = None
        energy_error_per_atom: float | None = None
        force_rmse: float | None = None
        force_mae: float | None = None
        force_max_error: float | None = None
        force_sse: float | None = None
        n_force_components = 0
        errors: list[str] = []

        if ref_energy is not None:
            try:
                pred_energy = float(pred_atoms.get_potential_energy())
                energy_error_per_atom = (pred_energy - ref_energy) / len(atoms)
            except Exception as exc:  # noqa: BLE001 - record frame failure
                errors.append(f"energy: {exc!r}")

        if ref_forces is not None:
            try:
                pred_forces = np.asarray(pred_atoms.get_forces(), dtype=float)
                if pred_forces.shape != ref_forces.shape:
                    raise ValueError(
                        f"{pred_forces.shape=} != reference {ref_forces.shape=}"
                    )
                force_errors = pred_forces - ref_forces
                abs_force_errors = np.abs(force_errors)
                force_sse = float(np.square(force_errors).sum())
                n_force_components = int(force_errors.size)
                force_rmse = float(np.sqrt(force_sse / n_force_components))
                force_mae = float(abs_force_errors.mean())
                force_max_error = float(abs_force_errors.max())
            except Exception as exc:  # noqa: BLE001 - record frame failure
                errors.append(f"forces: {exc!r}")

        rows.append(
            asdict(
                MdReferenceEvalResult(
                    model_key=model_key,
                    input_file=str(input_file),
                    system_name=system_name,
                    frame_idx=frame_idx,
                    n_atoms=len(atoms),
                    temperature_k=temperature_k,
                    ref_energy=ref_energy,
                    pred_energy=pred_energy,
                    energy_error_per_atom=energy_error_per_atom,
                    force_rmse=force_rmse,
                    force_mae=force_mae,
                    force_max_error=force_max_error,
                    force_sse=force_sse,
                    n_force_components=n_force_components,
                    status="failed" if errors else "completed",
                    error="; ".join(errors) or None,
                )
            )
        )
        if progress_bar is not None:
            progress_bar.update()

    return pd.DataFrame(rows)


def count_trajectory_frames(file_path: str | Path) -> int:
    """Count frames in an ASE-readable trajectory without keeping them all in memory."""
    return sum(1 for _frame in ase.io.iread(file_path, index=":"))


def calc_md_reference_metrics(df_eval: pd.DataFrame) -> dict[str, float]:
    """Aggregate reference-frame MD prediction errors for one model."""
    status = df_eval.get("status", pd.Series(dtype=str)).fillna("")

    def numeric_col(col_name: str) -> pd.Series:
        if col_name not in df_eval:
            return pd.Series(dtype=float)
        return pd.to_numeric(df_eval[col_name], errors="coerce")

    frame_idx = numeric_col("frame_idx")
    energy_errors = numeric_col("energy_error_per_atom").dropna()
    force_sse = numeric_col("force_sse")
    force_components = numeric_col("n_force_components").fillna(0)
    force_mask = force_sse.notna() & (force_components > 0)

    energy_rmse = (
        float(np.sqrt(np.mean(np.square(energy_errors))))
        if len(energy_errors) > 0
        else float("nan")
    )
    total_force_components = int(force_components[force_mask].sum())
    force_rmse = (
        float(np.sqrt(force_sse[force_mask].sum() / total_force_components))
        if total_force_components > 0
        else float("nan")
    )

    return {
        "energy_rmse": energy_rmse,
        "force_rmse": force_rmse,
        "n_frames": int(frame_idx.notna().sum()),
        "n_energy_frames": len(energy_errors),
        "n_force_frames": int(force_mask.sum()),
        "n_failed": int(status.eq("failed").sum()),
        "n_skipped": int(status.str.startswith("skipped").sum()),
    }


def _normalize_system_name(name: str) -> str:
    """Normalize system names for matching VDOS settings CSV entries."""
    return "".join(ch.lower() for ch in name if ch.isalnum())


def parse_system_temperature_from_dirname(dirname: str) -> tuple[str, int] | None:
    """Parse names like ``bulkAg_600K_Kapil`` into normalized system and temp."""
    match = re.match(r"^(?P<system>[^_]+)_(?P<temp>\d+)K(?:\b|[_-].*)$", dirname)
    if match is None:
        return None
    return _normalize_system_name(match["system"]), int(match["temp"])


def load_rdf_time_step_settings(
    path: str | Path | None,
) -> dict[tuple[str, int], dict[str, float | int]]:
    """Load per-system RDF/VDOS frame spacing settings from CSV."""
    if path is None:
        return {}

    settings: dict[tuple[str, int], dict[str, float | int]] = {}
    with open(path, newline="", encoding="utf-8") as file:
        for row in csv.DictReader(file):
            system = _normalize_system_name(row["System"])
            temperature = int(float(row["temperature"]))
            settings[(system, temperature)] = {
                "stride": int(float(row["stride"])),
                "dt_fs": float(row["dt"]),
                "padding": int(float(row["padding"])),
            }
    return settings


def _lcm_int(a: int, b: int) -> int:
    """Least common multiple for integers."""
    return abs(a * b) // math.gcd(a, b)


def _lcm_fraction(a: Fraction, b: Fraction) -> Fraction:
    """Least common multiple for positive rational numbers."""
    return Fraction(
        _lcm_int(a.numerator, b.numerator),
        math.gcd(a.denominator, b.denominator),
    )


def matched_frame_counts(
    *,
    n_ref_total: int,
    n_mlip_total: int,
    ref_dt_fs: float,
    mlip_dt_fs: float,
) -> tuple[int, int, float]:
    """Find frame counts spanning the same physical simulation time.

    Returns ``(n_ref_use, n_mlip_use, matched_time_fs)`` so that
    ``n_ref_use * ref_dt_fs == n_mlip_use * mlip_dt_fs`` and the matched time
    does not exceed either available trajectory length.
    """
    if n_ref_total <= 0 or n_mlip_total <= 0:
        return 0, 0, 0.0

    ref_dt = Fraction(str(ref_dt_fs))
    mlip_dt = Fraction(str(mlip_dt_fs))
    common_time = _lcm_fraction(ref_dt, mlip_dt)
    ref_time_max = n_ref_total * ref_dt
    mlip_time_max = n_mlip_total * mlip_dt
    n_common = min(ref_time_max // common_time, mlip_time_max // common_time)
    if n_common <= 0:
        return 0, 0, 0.0

    matched_time = n_common * common_time
    return int(matched_time / ref_dt), int(matched_time / mlip_dt), float(matched_time)


def _frame_dt_from_settings(
    system_name: str,
    settings: dict[tuple[str, int], dict[str, float | int]],
) -> float | None:
    """Look up frame spacing in fs for one system from RDF settings."""
    system_key = parse_system_temperature_from_dirname(system_name)
    if system_key is None:
        return None
    if config := settings.get(system_key):
        return float(config["dt_fs"])
    return None


def _vdos_parameters_from_settings(
    *,
    system_name: str,
    settings: dict[tuple[str, int], dict[str, float | int]],
    frame_dt_fs: float | None,
    stride: int | None,
    padding: int | None,
    label: str,
) -> tuple[int, float, int]:
    """Resolve trajectory stride, effective frame spacing, and FFT padding."""
    system_key = parse_system_temperature_from_dirname(system_name)
    system_settings = settings.get(system_key) if system_key is not None else None
    native_dt_fs = frame_dt_fs
    if native_dt_fs is None and system_settings is not None:
        native_dt_fs = float(system_settings["dt_fs"])
    if native_dt_fs is None:
        raise ValueError(
            f"Missing {label} VDOS frame spacing for {system_name}. "
            f"Pass {label}_frame_dt_fs or a matching settings CSV."
        )

    resolved_stride = stride
    if resolved_stride is None and system_settings is not None:
        resolved_stride = int(system_settings["stride"])
    resolved_stride = 1 if resolved_stride is None else resolved_stride

    resolved_padding = padding
    if resolved_padding is None and system_settings is not None:
        resolved_padding = int(system_settings["padding"])
    resolved_padding = 1 if resolved_padding is None else resolved_padding

    if native_dt_fs <= 0:
        raise ValueError(f"{label} VDOS frame spacing must be positive")
    if resolved_stride <= 0:
        raise ValueError(f"{label} VDOS stride must be positive")
    if resolved_padding <= 0:
        raise ValueError(f"{label} VDOS padding must be positive")

    return resolved_stride, native_dt_fs * resolved_stride, resolved_padding


def _atoms_to_mdtraj(frames: Sequence[Atoms]) -> Any:  # noqa: ANN401
    """Convert ASE frames to an MDTraj trajectory."""
    import mdtraj as mdt

    xyz = np.array([atoms.get_positions() for atoms in frames], dtype=np.float64) / 10
    symbols = frames[0].get_chemical_symbols()
    topology = mdt.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("SYS", chain)
    for symbol in symbols:
        topology.add_atom(
            symbol,
            element=mdt.element.get_by_symbol(symbol),
            residue=residue,
        )

    md_traj = mdt.Trajectory(xyz=xyz, topology=topology)
    cells = (
        np.array([atoms.get_cell().array for atoms in frames], dtype=np.float64) / 10
    )
    md_traj.unitcell_vectors = cells
    md_traj.unitcell_lengths = np.linalg.norm(cells, axis=2)
    return md_traj


def _compute_rdf_with_mdtraj(
    frames: Sequence[Atoms],
    *,
    nbins: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute RDF with MDTraj, matching the source RDF scripts."""
    md_traj = _atoms_to_mdtraj(frames)
    atoms = np.arange(md_traj.n_atoms)
    pairs = md_traj.top.select_pairs(atoms, atoms)
    cell_lengths_ang = md_traj.unitcell_lengths[0] * 10
    r_max_nm = np.min(cell_lengths_ang) / 20

    import mdtraj as mdt

    r_nm, g_r = mdt.compute_rdf(
        md_traj,
        pairs=pairs,
        r_range=(0, r_max_nm),
        n_bins=nbins,
        periodic=True,
    )
    return r_nm * 10, g_r


def _compute_rdf_with_ase(
    frames: Sequence[Atoms],
    *,
    nbins: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute an all-pair, time-averaged RDF using ASE minimum-image distances."""
    if len(frames) == 0:
        raise ValueError("Cannot compute RDF for an empty trajectory")

    n_atoms = len(frames[0])
    if n_atoms < 2:
        raise ValueError("Cannot compute RDF with fewer than two atoms")

    cell_lengths_ang = frames[0].cell.lengths()
    r_max = float(np.min(cell_lengths_ang) / 2)
    bin_edges = np.linspace(0, r_max, nbins + 1)
    counts = np.zeros(nbins, dtype=float)
    upper_tri = np.triu_indices(n_atoms, k=1)
    volumes: list[float] = []

    for atoms in frames:
        if len(atoms) != n_atoms:
            raise ValueError("All RDF trajectory frames must have the same atom count")
        distances = atoms.get_all_distances(mic=True)[upper_tri]
        counts += np.histogram(distances, bins=bin_edges)[0]
        volumes.append(float(atoms.get_volume()))

    r = (bin_edges[:-1] + bin_edges[1:]) / 2
    shell_volumes = 4 / 3 * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)
    density = n_atoms / float(np.mean(volumes))
    normalizer = len(frames) * n_atoms * density * shell_volumes / 2
    g_r = np.divide(
        counts,
        normalizer,
        out=np.zeros_like(counts),
        where=normalizer > 0,
    )
    return r, g_r


def compute_rdf(
    frames: Sequence[Atoms],
    *,
    nbins: int = 500,
    prefer_mdtraj: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute time-averaged all-pair RDF for ASE trajectory frames.

    If MDTraj is installed, this follows the implementation used by the source
    benchmark scripts. A pure ASE fallback keeps tests and lightweight installs usable.
    """
    if prefer_mdtraj:
        try:
            return _compute_rdf_with_mdtraj(frames, nbins=nbins)
        except ModuleNotFoundError:
            pass
    return _compute_rdf_with_ase(frames, nbins=nbins)


def rdf_error(
    ref_rdf: tuple[np.ndarray, np.ndarray],
    test_rdf: tuple[np.ndarray, np.ndarray],
) -> float:
    """Return source-compatible RDF error in percent."""
    r_ref, g_ref = ref_rdf
    r_test, g_test = test_rdf
    g_test_interp = np.interp(r_ref, r_test, g_test)
    denominator = float(np.sum(np.abs(g_ref - 1)))
    if denominator == 0:
        return 100.0
    numerator = float(np.sum(np.abs(g_ref - g_test_interp)))
    return min(1.0, numerator / denominator) * 100


def save_rdf_csv(r: np.ndarray, g_r: np.ndarray, path: str | Path) -> None:
    """Save RDF curve as ``r_A,g_r`` CSV."""
    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        file_path,
        np.column_stack([r, g_r]),
        delimiter=",",
        header="r_A,g_r",
        comments="",
    )


def compute_vdos_hann(
    frames: Sequence[Atoms],
    *,
    frame_dt_fs: float,
    padding: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute VDOS using the Hann/FFT algorithm used by the final scripts.

    This mirrors ``get_VDOS_padding_hann_og.py`` from the source workflow:
    velocities are finite differences of Cartesian positions, an independently
    normalized ACF and power spectrum is computed for each degree of freedom,
    and the spectra are summed.
    """
    from scipy import signal

    if len(frames) < 2:
        raise ValueError("At least two trajectory frames are required for VDOS")
    if frame_dt_fs <= 0:
        raise ValueError("VDOS frame spacing must be positive")
    if padding <= 0:
        raise ValueError("VDOS padding must be positive")

    n_steps = len(frames)
    n_atoms = len(frames[0])
    coordinates = np.empty((n_steps, n_atoms, 3), dtype=float)
    for idx, atoms in enumerate(frames):
        if len(atoms) != n_atoms:
            raise ValueError("All VDOS trajectory frames must have the same atom count")
        coordinates[idx] = atoms.get_positions()

    delta_t_s = frame_dt_fs * 1e-15
    window = signal.windows.hann(n_steps, sym=False)
    window_efficiency = float(np.sum(window) / n_steps)
    if window_efficiency == 0:
        raise ValueError("VDOS Hann window has zero efficiency")
    window = window / window_efficiency
    n_fft = padding * int(2 ** math.ceil(math.log2(n_steps)))
    total_intensity = np.zeros(n_fft, dtype=float)

    for atom_idx in range(n_atoms):
        for component_idx in range(3):
            velocities = (
                np.gradient(coordinates[:, atom_idx, component_idx]) / delta_t_s
            )
            centered = velocities - np.mean(velocities)
            norm = float(np.sum(np.square(centered)))
            if norm == 0:
                acf = np.zeros_like(centered)
            else:
                acf = (
                    signal.fftconvolve(
                        centered,
                        centered[::-1],
                        mode="full",
                    )[n_steps - 1 :]
                    / norm
                )
            component_fft = np.fft.fft(acf * window, n=n_fft) / n_steps
            total_intensity += np.square(np.abs(component_fft))

    half = n_fft // 2
    wavenumber_cm = np.fft.fftfreq(n_fft, delta_t_s * SPEED_OF_LIGHT_CM_S)[:half]
    return wavenumber_cm, total_intensity[:half]


def save_vdos_dat(
    wavenumber_cm: np.ndarray,
    intensity: np.ndarray,
    path: str | Path,
) -> None:
    """Save a VDOS spectrum in the two-column format of the source scripts."""
    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        file_path,
        np.column_stack([wavenumber_cm, intensity]),
        header="# Wavenumber(cm-1)   Intensity(a.u.)",
        comments="",
    )


def _vdos_spectrum_cm_to_ev(
    spectrum: tuple[np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Clean and convert a VDOS spectrum wavenumber axis from cm^-1 to eV."""
    wavenumber_cm, intensity = spectrum
    wavenumber_cm = np.asarray(wavenumber_cm, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    valid = np.isfinite(wavenumber_cm) & np.isfinite(intensity) & (wavenumber_cm >= 0)
    energy_ev = wavenumber_cm[valid] * CM_INV_TO_EV
    intensity = intensity[valid]
    if energy_ev.size < 2:
        raise ValueError("VDOS spectrum must contain at least two valid points")

    order = np.argsort(energy_ev)
    energy_ev = energy_ev[order]
    intensity = intensity[order]
    unique_energy, unique_indices = np.unique(energy_ev, return_index=True)
    return unique_energy, intensity[unique_indices]


def _trapz_integral(y: np.ndarray, x: np.ndarray) -> float:
    """Integrate an array using the NumPy version available at runtime."""
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x=x))
    legacy_trapz: Any = vars(np)["trapz"]
    return float(legacy_trapz(y, x=x))


def vdos_error(
    ref_vdos: tuple[np.ndarray, np.ndarray],
    test_vdos: tuple[np.ndarray, np.ndarray],
    *,
    e_min: float | None = None,
    e_max: float | None = None,
    normalize_area: bool = True,
) -> float:
    """Return normalized VDOS error on the source workflow's 0..1 scale."""
    ref_energy, ref_intensity = _vdos_spectrum_cm_to_ev(ref_vdos)
    test_energy, test_intensity = _vdos_spectrum_cm_to_ev(test_vdos)
    lower = max(float(ref_energy.min()), float(test_energy.min()))
    upper = min(float(ref_energy.max()), float(test_energy.max()))
    if e_min is not None:
        lower = max(lower, e_min)
    if e_max is not None:
        upper = min(upper, e_max)
    if not np.isfinite(lower) or not np.isfinite(upper) or upper <= lower:
        return float("nan")

    grid = np.unique(
        np.concatenate(
            [
                ref_energy[(ref_energy >= lower) & (ref_energy <= upper)],
                test_energy[(test_energy >= lower) & (test_energy <= upper)],
            ]
        )
    )
    if grid.size < 2:
        return float("nan")

    ref_interp = np.clip(np.interp(grid, ref_energy, ref_intensity), 0, None)
    test_interp = np.clip(np.interp(grid, test_energy, test_intensity), 0, None)
    if normalize_area:
        ref_area = _trapz_integral(ref_interp, grid)
        test_area = _trapz_integral(test_interp, grid)
        if ref_area <= 0 or test_area <= 0:
            raise ValueError("VDOS spectrum has non-positive area; cannot normalize")
        ref_interp = ref_interp / ref_area
        test_interp = test_interp / test_area

    numerator = _trapz_integral(np.abs(ref_interp - test_interp), grid)
    denominator = _trapz_integral(ref_interp, grid) + _trapz_integral(test_interp, grid)
    if denominator <= 0 or not np.isfinite(denominator):
        return float("nan")
    return numerator / denominator


def _vdos_curve_paths(
    *,
    output_dir: str | Path,
    mode: VdosMode,
    model_key: str,
    system_name: str,
    ref_stride: int,
    mlip_stride: int,
) -> tuple[Path, Path]:
    """Return VDOS curve paths matching the source result layout."""
    mode_name = (
        "vdos_results_hann_same_simulation_length"
        if mode == "same"
        else "vdos_results_hann_different_simulation_length"
    )
    root = Path(output_dir) / mode_name
    ref_stem = "traj" if ref_stride == 1 else f"traj_stride{ref_stride}"
    mlip_stem = f"nvt_{model_key}"
    if mlip_stride > 1:
        mlip_stem = f"{mlip_stem}_stride{mlip_stride}"

    if mode == "same":
        ref_path = (
            root
            / "reference_matched"
            / system_name
            / f"{ref_stem}_match_nvt_{model_key}_vdos_hann.dat"
        )
    else:
        ref_path = root / "reference" / system_name / f"{ref_stem}_vdos_hann.dat"
    mlip_path = root / "mlip" / system_name / f"{mlip_stem}_vdos_hann.dat"
    return ref_path, mlip_path


def _stress_as_3x3(stress: np.ndarray) -> np.ndarray:
    """Normalize an ASE stress tensor into a 3 by 3 matrix."""
    array = np.asarray(stress, dtype=float)
    if array.shape == (3, 3):
        return array
    flat = array.reshape(-1)
    if flat.size == 6:
        xx, yy, zz, yz, xz, xy = flat.tolist()
        return np.array(
            [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]],
            dtype=float,
        )
    if flat.size == 9:
        return flat.reshape(3, 3)
    raise ValueError(f"Unsupported stress shape: {array.shape}")


def pressure_from_stress_gpa(
    stress: np.ndarray,
    *,
    pressure_mode: Literal["3d", "2d"] = "3d",
    plane_2d: Literal["xy", "xz", "yz"] = "xy",
) -> float:
    """Convert ASE stress in eV/Angstrom^3 into pressure in GPa."""
    stress_matrix = _stress_as_3x3(stress)
    if pressure_mode == "3d":
        sigma = float(np.trace(stress_matrix) / 3)
    elif plane_2d == "xy":
        sigma = float((stress_matrix[0, 0] + stress_matrix[1, 1]) / 2)
    elif plane_2d == "xz":
        sigma = float((stress_matrix[0, 0] + stress_matrix[2, 2]) / 2)
    else:
        sigma = float((stress_matrix[1, 1] + stress_matrix[2, 2]) / 2)
    return -sigma * EV_PER_A3_TO_GPA


def _stored_frame_stress(atoms: Atoms) -> np.ndarray:
    """Read stress stored on a reference trajectory frame."""
    try:
        return np.asarray(atoms.get_stress(voigt=False), dtype=float)
    except Exception:  # noqa: BLE001 - fall back to stored stress metadata
        if "stress" in atoms.info:
            return np.asarray(atoms.info["stress"], dtype=float)
    raise KeyError("No stress found in reference frame")


def pressure_histogram_error(
    ref_values: np.ndarray,
    mlip_values: np.ndarray,
    *,
    bins: int = 80,
) -> dict[str, float]:
    """Compute the source pressure histogram percentage error."""
    if bins < 2:
        raise ValueError("Pressure histogram bins must be at least 2")
    ref_values = np.asarray(ref_values, dtype=float)
    mlip_values = np.asarray(mlip_values, dtype=float)
    if ref_values.size == 0 or mlip_values.size == 0:
        raise ValueError("Cannot compute pressure histogram error for empty values")

    lower = float(min(np.min(ref_values), np.min(mlip_values)))
    upper = float(max(np.max(ref_values), np.max(mlip_values)))
    if not np.isfinite(lower) or not np.isfinite(upper):
        lower, upper = -1.0, 1.0
    if upper <= lower:
        upper = lower + 1e-6
    edges = np.linspace(lower, upper, bins + 1)
    widths = np.diff(edges)

    ref_hist, _ = np.histogram(ref_values, bins=edges, density=True)
    mlip_hist, _ = np.histogram(mlip_values, bins=edges, density=True)
    ref_area = float(np.sum(ref_hist * widths))
    mlip_area = float(np.sum(mlip_hist * widths))
    if ref_area <= 0 or mlip_area <= 0:
        raise ValueError("Pressure histogram has non-positive area")

    ref_hist = ref_hist / ref_area
    mlip_hist = mlip_hist / mlip_area
    ref_area = float(np.sum(ref_hist * widths))
    mlip_area = float(np.sum(mlip_hist * widths))
    numerator = float(np.sum(np.abs(ref_hist - mlip_hist) * widths))
    denominator = ref_area + mlip_area
    distance = numerator / denominator
    error_fraction = float(np.clip(distance, 0, 1))
    return {
        "pressure_error_percent": 100 * error_fraction,
        "pressure_histogram_l1_area": numerator,
        "pressure_histogram_distance": distance,
        "pressure_ref_histogram_area": ref_area,
        "pressure_mlip_histogram_area": mlip_area,
        "pressure_histogram_min_GPa": lower,
        "pressure_histogram_max_GPa": upper,
    }


def _pressure_value_paths(
    *,
    output_dir: str | Path,
    model_key: str,
    system_name: str,
) -> tuple[Path, Path]:
    """Return per-frame pressure value paths for one matched trajectory pair."""
    root = Path(output_dir) / "pressure_results_same_simulation_length"
    return (
        root / "reference" / f"{system_name}__{model_key}_pressure_per_frame.csv",
        root / "mlip" / f"{system_name}__{model_key}_pressure_per_frame.csv",
    )


def _save_pressure_values(
    values: np.ndarray,
    *,
    path: Path,
    system_name: str,
    trajectory_file: str | Path,
    column_name: str,
    pressure_mode: str,
    plane_used: str | None,
) -> None:
    """Write per-frame pressure values in a source-compatible CSV shape."""
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "system": system_name,
            "trajectory_file": str(trajectory_file),
            "frame_index": np.arange(len(values)),
            column_name: values,
            "pressure_mode": pressure_mode,
            "plane_used": plane_used,
        }
    ).to_csv(path, index=False)


def _rdf_curve_paths(
    *,
    output_dir: str | Path,
    mode: RdfMode,
    model_key: str,
    system_name: str,
) -> tuple[Path, Path]:
    """Return reference and MLIP RDF curve paths matching source script layout."""
    mode_name = (
        "rdf_same_simulation_length_saved"
        if mode == "same"
        else "rdf_different_simulation_length_saved"
    )
    root = Path(output_dir) / mode_name
    if mode == "same":
        ref_name = f"{system_name}__{model_key}.csv"
        mlip_name = f"{system_name}__{model_key}.csv"
    else:
        ref_name = f"{system_name}.csv"
        mlip_name = f"{system_name}.csv"
    return root / "reference" / ref_name, root / "mlip" / model_key / mlip_name


def find_mlip_trajectory(
    *,
    mlip_dir: str | Path,
    system_name: str,
    model_key: str,
) -> Path | None:
    """Find an MLIP-predicted trajectory for a reference system."""
    root = Path(mlip_dir)
    if root.is_file():
        return root

    candidates = (
        root / system_name / model_key / "simulation.pos_0.extxyz",
        root / system_name / model_key / "simulation.pos_0.xyz",
        root / system_name / f"nvt_{model_key}.extxyz",
        root / system_name / f"nvt_{model_key}.xyz",
        root / system_name / "traj.extxyz",
        root / system_name / "traj.xyz",
        root / f"nvt_{model_key}.extxyz",
        root / f"nvt_{model_key}.xyz",
        root / f"{system_name}.extxyz",
        root / f"{system_name}.xyz",
    )
    return next((path for path in candidates if path.is_file()), None)


def evaluate_rdf_trajectory_pair(
    *,
    model_key: str,
    system_name: str,
    ref_file: str | Path,
    mlip_file: str | Path,
    output_dir: str | Path,
    config: MdRdfConfig,
    ref_frame_dt_fs: float | None = None,
    mlip_frame_dt_fs: float | None = None,
) -> MdRdfEvalResult:
    """Evaluate RDF error for one reference/MLIP trajectory pair."""
    ref_traj_all = read_trajectory(ref_file)
    mlip_traj_all = read_trajectory(mlip_file)
    if len(ref_traj_all) == 0:
        raise ValueError(f"Empty reference trajectory: {ref_file}")
    if len(mlip_traj_all) == 0:
        raise ValueError(f"Empty MLIP trajectory: {mlip_file}")

    matched_time_fs: float | None = None
    ref_traj = ref_traj_all
    mlip_traj = mlip_traj_all
    if config.mode == "same":
        if ref_frame_dt_fs is None or mlip_frame_dt_fs is None:
            n_use = min(len(ref_traj_all), len(mlip_traj_all))
            ref_traj = ref_traj_all[:n_use]
            mlip_traj = mlip_traj_all[:n_use]
        else:
            n_ref_use, n_mlip_use, matched_time_fs = matched_frame_counts(
                n_ref_total=len(ref_traj_all),
                n_mlip_total=len(mlip_traj_all),
                ref_dt_fs=ref_frame_dt_fs,
                mlip_dt_fs=mlip_frame_dt_fs,
            )
            if n_ref_use <= 0 or n_mlip_use <= 0:
                raise ValueError("Could not find positive matched simulation time")
            ref_traj = ref_traj_all[:n_ref_use]
            mlip_traj = mlip_traj_all[:n_mlip_use]

    ref_rdf = compute_rdf(ref_traj, nbins=config.nbins)
    mlip_rdf = compute_rdf(mlip_traj, nbins=config.nbins)
    error = rdf_error(ref_rdf, mlip_rdf)

    if config.write_rdf_curves:
        ref_rdf_path, mlip_rdf_path = _rdf_curve_paths(
            output_dir=output_dir,
            mode=config.mode,
            model_key=model_key,
            system_name=system_name,
        )
        save_rdf_csv(ref_rdf[0], ref_rdf[1], ref_rdf_path)
        save_rdf_csv(mlip_rdf[0], mlip_rdf[1], mlip_rdf_path)

    return MdRdfEvalResult(
        model_key=model_key,
        system_name=system_name,
        ref_file=str(ref_file),
        mlip_file=str(mlip_file),
        rdf_error=error,
        n_ref_frames=len(ref_traj),
        n_mlip_frames=len(mlip_traj),
        ref_frame_dt_fs=ref_frame_dt_fs,
        mlip_frame_dt_fs=mlip_frame_dt_fs,
        matched_time_fs=matched_time_fs,
    )


def calc_md_rdf_metrics(df_eval: pd.DataFrame) -> dict[str, float]:
    """Aggregate source-compatible RDF error rows for one model."""
    status = df_eval.get("status", pd.Series(dtype=str)).fillna("")
    if "rdf_error" in df_eval:
        errors = pd.to_numeric(df_eval["rdf_error"], errors="coerce")
    elif "RDF_Error" in df_eval:
        errors = pd.to_numeric(df_eval["RDF_Error"], errors="coerce")
    else:
        raise ValueError("RDF metrics require an rdf_error or RDF_Error column")
    errors = errors.dropna()
    return {
        "rdf_error": float(errors.mean()) if len(errors) else float("nan"),
        "n_rdf_systems": len(errors),
        "n_rdf_failed": int(status.eq("failed").sum()),
        "n_rdf_skipped": int(status.str.startswith("skipped").sum()),
    }


def evaluate_rdf_dataset(
    *,
    model_key: str,
    ref_dir: str | Path,
    mlip_dir: str | Path,
    output_dir: str | Path,
    config: MdRdfConfig | None = None,
    md_config: MdRunConfig | None = None,
    extensions: Sequence[str] = (".extxyz", ".xyz", ".extxyz.gz", ".xyz.gz"),
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Evaluate RDF error for reference trajectories and MLIP rollouts."""
    config = config or MdRdfConfig()
    md_config = md_config or MdRunConfig()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    ref_settings = load_rdf_time_step_settings(config.ref_settings_path)
    mlip_settings = load_rdf_time_step_settings(config.mlip_settings_path)
    rows: list[dict[str, Any]] = []
    all_ref_files = find_trajectory_files(ref_dir, extensions=extensions)

    with tqdm(
        total=len(all_ref_files),
        desc=f"Evaluating RDF with {model_key}",
        unit="system",
    ) as progress_bar:
        for ref_file in all_ref_files:
            system_name = ref_file.parent.name
            progress_bar.set_postfix_str(system_name, refresh=False)

            if should_skip_trajectory(ref_file, md_config):
                rows.append(
                    asdict(
                        MdRdfEvalResult(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=str(ref_file),
                            mlip_file=None,
                            rdf_error=None,
                            n_ref_frames=0,
                            n_mlip_frames=0,
                            ref_frame_dt_fs=None,
                            mlip_frame_dt_fs=None,
                            matched_time_fs=None,
                            status="skipped_input",
                        )
                    )
                )
                progress_bar.update()
                continue

            mlip_file = find_mlip_trajectory(
                mlip_dir=mlip_dir,
                system_name=system_name,
                model_key=model_key,
            )
            if mlip_file is None:
                rows.append(
                    asdict(
                        MdRdfEvalResult(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=str(ref_file),
                            mlip_file=None,
                            rdf_error=None,
                            n_ref_frames=0,
                            n_mlip_frames=0,
                            ref_frame_dt_fs=None,
                            mlip_frame_dt_fs=None,
                            matched_time_fs=None,
                            status="skipped_missing_mlip",
                            error=f"No MLIP trajectory found in {mlip_dir}",
                        )
                    )
                )
                progress_bar.update()
                continue

            ref_frame_dt_fs = config.ref_frame_dt_fs or _frame_dt_from_settings(
                system_name, ref_settings
            )
            mlip_frame_dt_fs = config.mlip_frame_dt_fs or _frame_dt_from_settings(
                system_name, mlip_settings
            )

            try:
                rows.append(
                    asdict(
                        evaluate_rdf_trajectory_pair(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=ref_file,
                            mlip_file=mlip_file,
                            output_dir=output_path,
                            config=config,
                            ref_frame_dt_fs=ref_frame_dt_fs,
                            mlip_frame_dt_fs=mlip_frame_dt_fs,
                        )
                    )
                )
            except Exception as exc:  # noqa: BLE001 - record per-system failure
                rows.append(
                    asdict(
                        MdRdfEvalResult(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=str(ref_file),
                            mlip_file=str(mlip_file),
                            rdf_error=None,
                            n_ref_frames=0,
                            n_mlip_frames=0,
                            ref_frame_dt_fs=ref_frame_dt_fs,
                            mlip_frame_dt_fs=mlip_frame_dt_fs,
                            matched_time_fs=None,
                            status="failed",
                            error=f"{exc!r}\n{traceback.format_exc()}",
                        )
                    )
                )
            progress_bar.update()

    df_eval = pd.DataFrame(rows)
    mode_name = "same_simulation_length" if config.mode == "same" else "different"
    pred_path = output_path / f"md_rdf_{mode_name}_{model_key}.csv"
    df_eval.to_csv(pred_path, index=False)

    metrics = calc_md_rdf_metrics(df_eval)
    with open(
        output_path / f"md_rdf_metrics_{mode_name}_{model_key}.json",
        mode="w",
        encoding="utf-8",
    ) as file:
        json.dump(metrics, file, indent=2)

    return df_eval, metrics


def evaluate_vdos_trajectory_pair(
    *,
    model_key: str,
    system_name: str,
    ref_file: str | Path,
    mlip_file: str | Path,
    output_dir: str | Path,
    config: MdVdosConfig,
    ref_frame_dt_fs: float,
    mlip_frame_dt_fs: float,
    ref_stride: int,
    mlip_stride: int,
    ref_padding: int,
    mlip_padding: int,
) -> MdVdosEvalResult:
    """Evaluate source-compatible VDOS error for one trajectory pair."""
    ref_traj = read_trajectory(ref_file)[::ref_stride]
    mlip_traj = read_trajectory(mlip_file)[::mlip_stride]
    if not ref_traj:
        raise ValueError(f"Empty reference trajectory after stride: {ref_file}")
    if not mlip_traj:
        raise ValueError(f"Empty MLIP trajectory after stride: {mlip_file}")

    matched_time_fs: float | None = None
    if config.mode == "same":
        if not math.isclose(
            ref_frame_dt_fs,
            mlip_frame_dt_fs,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError(
                "Same-length VDOS requires equal effective frame spacing after "
                f"stride: reference={ref_frame_dt_fs} fs, "
                f"MLIP={mlip_frame_dt_fs} fs"
            )
        n_use = min(len(ref_traj), len(mlip_traj))
        ref_traj = ref_traj[:n_use]
        mlip_traj = mlip_traj[:n_use]
        matched_time_fs = n_use * ref_frame_dt_fs

    ref_vdos = compute_vdos_hann(
        ref_traj,
        frame_dt_fs=ref_frame_dt_fs,
        padding=ref_padding,
    )
    mlip_vdos = compute_vdos_hann(
        mlip_traj,
        frame_dt_fs=mlip_frame_dt_fs,
        padding=mlip_padding,
    )
    error = vdos_error(
        ref_vdos,
        mlip_vdos,
        e_min=config.e_min,
        e_max=config.e_max,
        normalize_area=config.normalize_area,
    )

    if config.write_vdos_curves:
        ref_vdos_path, mlip_vdos_path = _vdos_curve_paths(
            output_dir=output_dir,
            mode=config.mode,
            model_key=model_key,
            system_name=system_name,
            ref_stride=ref_stride,
            mlip_stride=mlip_stride,
        )
        save_vdos_dat(*ref_vdos, ref_vdos_path)
        save_vdos_dat(*mlip_vdos, mlip_vdos_path)

    return MdVdosEvalResult(
        model_key=model_key,
        system_name=system_name,
        ref_file=str(ref_file),
        mlip_file=str(mlip_file),
        vdos_error=error,
        n_ref_frames=len(ref_traj),
        n_mlip_frames=len(mlip_traj),
        ref_frame_dt_fs=ref_frame_dt_fs,
        mlip_frame_dt_fs=mlip_frame_dt_fs,
        ref_stride=ref_stride,
        mlip_stride=mlip_stride,
        ref_padding=ref_padding,
        mlip_padding=mlip_padding,
        matched_time_fs=matched_time_fs,
    )


def calc_md_vdos_metrics(df_eval: pd.DataFrame) -> dict[str, float]:
    """Aggregate normalized VDOS error rows for one model."""
    status = df_eval.get("status", pd.Series(dtype=str)).fillna("")
    if "vdos_error" in df_eval:
        errors = pd.to_numeric(df_eval["vdos_error"], errors="coerce")
    elif "final_mean_vdos_error" in df_eval:
        errors = pd.to_numeric(df_eval["final_mean_vdos_error"], errors="coerce")
    elif "mean_vdos_error" in df_eval:
        errors = pd.to_numeric(df_eval["mean_vdos_error"], errors="coerce")
    else:
        raise ValueError(
            "VDOS metrics require a vdos_error, final_mean_vdos_error, or "
            "mean_vdos_error column"
        )
    errors = errors.dropna()
    return {
        "vdos_error": float(errors.mean()) if len(errors) else float("nan"),
        "n_vdos_systems": len(errors),
        "n_vdos_failed": int(status.eq("failed").sum()),
        "n_vdos_skipped": int(status.str.startswith("skipped").sum()),
    }


def evaluate_vdos_dataset(
    *,
    model_key: str,
    ref_dir: str | Path,
    mlip_dir: str | Path,
    output_dir: str | Path,
    config: MdVdosConfig | None = None,
    md_config: MdRunConfig | None = None,
    extensions: Sequence[str] = (".extxyz", ".xyz", ".extxyz.gz", ".xyz.gz"),
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Evaluate VDOS error for reference trajectories and MLIP rollouts."""
    config = config or MdVdosConfig()
    md_config = md_config or MdRunConfig()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    ref_settings = load_rdf_time_step_settings(config.ref_settings_path)
    mlip_settings = load_rdf_time_step_settings(config.mlip_settings_path)
    rows: list[dict[str, Any]] = []
    all_ref_files = find_trajectory_files(ref_dir, extensions=extensions)

    with tqdm(
        total=len(all_ref_files),
        desc=f"Evaluating VDOS with {model_key}",
        unit="system",
    ) as progress_bar:
        for ref_file in all_ref_files:
            system_name = ref_file.parent.name
            progress_bar.set_postfix_str(system_name, refresh=False)

            if should_skip_trajectory(ref_file, md_config):
                rows.append(
                    asdict(
                        MdVdosEvalResult(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=str(ref_file),
                            mlip_file=None,
                            vdos_error=None,
                            n_ref_frames=0,
                            n_mlip_frames=0,
                            ref_frame_dt_fs=None,
                            mlip_frame_dt_fs=None,
                            ref_stride=None,
                            mlip_stride=None,
                            ref_padding=None,
                            mlip_padding=None,
                            matched_time_fs=None,
                            status="skipped_input",
                        )
                    )
                )
                progress_bar.update()
                continue

            mlip_file = find_mlip_trajectory(
                mlip_dir=mlip_dir,
                system_name=system_name,
                model_key=model_key,
            )
            if mlip_file is None:
                rows.append(
                    asdict(
                        MdVdosEvalResult(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=str(ref_file),
                            mlip_file=None,
                            vdos_error=None,
                            n_ref_frames=0,
                            n_mlip_frames=0,
                            ref_frame_dt_fs=None,
                            mlip_frame_dt_fs=None,
                            ref_stride=None,
                            mlip_stride=None,
                            ref_padding=None,
                            mlip_padding=None,
                            matched_time_fs=None,
                            status="skipped_missing_mlip",
                            error=f"No MLIP trajectory found in {mlip_dir}",
                        )
                    )
                )
                progress_bar.update()
                continue

            try:
                ref_stride, ref_dt_fs, ref_padding = _vdos_parameters_from_settings(
                    system_name=system_name,
                    settings=ref_settings,
                    frame_dt_fs=config.ref_frame_dt_fs,
                    stride=config.ref_stride,
                    padding=config.ref_padding,
                    label="ref",
                )
                mlip_stride, mlip_dt_fs, mlip_padding = _vdos_parameters_from_settings(
                    system_name=system_name,
                    settings=mlip_settings,
                    frame_dt_fs=config.mlip_frame_dt_fs,
                    stride=config.mlip_stride,
                    padding=config.mlip_padding,
                    label="mlip",
                )
                rows.append(
                    asdict(
                        evaluate_vdos_trajectory_pair(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=ref_file,
                            mlip_file=mlip_file,
                            output_dir=output_path,
                            config=config,
                            ref_frame_dt_fs=ref_dt_fs,
                            mlip_frame_dt_fs=mlip_dt_fs,
                            ref_stride=ref_stride,
                            mlip_stride=mlip_stride,
                            ref_padding=ref_padding,
                            mlip_padding=mlip_padding,
                        )
                    )
                )
            except Exception as exc:  # noqa: BLE001 - record per-system failure
                rows.append(
                    asdict(
                        MdVdosEvalResult(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=str(ref_file),
                            mlip_file=str(mlip_file),
                            vdos_error=None,
                            n_ref_frames=0,
                            n_mlip_frames=0,
                            ref_frame_dt_fs=None,
                            mlip_frame_dt_fs=None,
                            ref_stride=None,
                            mlip_stride=None,
                            ref_padding=None,
                            mlip_padding=None,
                            matched_time_fs=None,
                            status="failed",
                            error=f"{exc!r}\n{traceback.format_exc()}",
                        )
                    )
                )
            progress_bar.update()

    df_eval = pd.DataFrame(rows)
    mode_name = (
        "same_simulation_length"
        if config.mode == "same"
        else "different_simulation_length"
    )
    pred_path = output_path / f"md_vdos_{mode_name}_{model_key}.csv"
    df_eval.to_csv(pred_path, index=False)

    metrics = calc_md_vdos_metrics(df_eval)
    with open(
        output_path / f"md_vdos_metrics_{mode_name}_{model_key}.json",
        mode="w",
        encoding="utf-8",
    ) as file:
        json.dump(metrics, file, indent=2)

    return df_eval, metrics


def evaluate_pressure_trajectory_pair(
    *,
    model_key: str,
    system_name: str,
    ref_file: str | Path,
    mlip_file: str | Path,
    output_dir: str | Path,
    calculator: Calculator,
    config: MdPressureConfig,
    ref_frame_dt_fs: float,
    mlip_frame_dt_fs: float,
) -> MdPressureEvalResult:
    """Evaluate matched-length pressure histogram error for one trajectory pair."""
    ref_traj = read_trajectory(ref_file)
    mlip_traj = read_trajectory(mlip_file)
    n_ref_use, n_mlip_use, matched_time_fs = matched_frame_counts(
        n_ref_total=len(ref_traj),
        n_mlip_total=len(mlip_traj),
        ref_dt_fs=ref_frame_dt_fs,
        mlip_dt_fs=mlip_frame_dt_fs,
    )
    if n_ref_use <= 0 or n_mlip_use <= 0:
        raise ValueError("Could not find positive matched simulation time for pressure")

    normalized_system = _normalize_system_name(system_name.split("_", maxsplit=1)[0])
    normalized_2d = {_normalize_system_name(name) for name in config.two_d_systems}
    is_2d = normalized_system in normalized_2d
    pressure_mode: Literal["3d", "2d"] = "2d" if is_2d else "3d"
    plane_used = config.plane_2d if is_2d else None

    ref_values = np.asarray(
        [
            pressure_from_stress_gpa(
                _stored_frame_stress(atoms),
                pressure_mode=pressure_mode,
                plane_2d=config.plane_2d,
            )
            for atoms in ref_traj[:n_ref_use]
        ]
    )
    mlip_values: list[float] = []
    for atoms in mlip_traj[:n_mlip_use]:
        predicted = atoms.copy()
        predicted.calc = calculator
        mlip_values.append(
            pressure_from_stress_gpa(
                predicted.get_stress(voigt=False),
                pressure_mode=pressure_mode,
                plane_2d=config.plane_2d,
            )
        )
    mlip_array = np.asarray(mlip_values)
    scores = pressure_histogram_error(ref_values, mlip_array, bins=config.bins)

    if config.write_pressure_values:
        ref_pressure_path, mlip_pressure_path = _pressure_value_paths(
            output_dir=output_dir,
            model_key=model_key,
            system_name=system_name,
        )
        _save_pressure_values(
            ref_values,
            path=ref_pressure_path,
            system_name=system_name,
            trajectory_file=ref_file,
            column_name="pressure_ref_GPa",
            pressure_mode=pressure_mode,
            plane_used=plane_used,
        )
        _save_pressure_values(
            mlip_array,
            path=mlip_pressure_path,
            system_name=system_name,
            trajectory_file=mlip_file,
            column_name="pressure_GPa",
            pressure_mode=pressure_mode,
            plane_used=plane_used,
        )

    return MdPressureEvalResult(
        model_key=model_key,
        system_name=system_name,
        ref_file=str(ref_file),
        mlip_file=str(mlip_file),
        bins=config.bins,
        n_ref_frames=len(ref_values),
        n_mlip_frames=len(mlip_array),
        ref_frame_dt_fs=ref_frame_dt_fs,
        mlip_frame_dt_fs=mlip_frame_dt_fs,
        matched_time_fs=matched_time_fs,
        pressure_mode=pressure_mode,
        plane_used=plane_used,
        pressure_error_percent=scores["pressure_error_percent"],
        pressure_histogram_l1_area=scores["pressure_histogram_l1_area"],
        pressure_histogram_distance=scores["pressure_histogram_distance"],
        pressure_ref_histogram_area=scores["pressure_ref_histogram_area"],
        pressure_mlip_histogram_area=scores["pressure_mlip_histogram_area"],
        pressure_histogram_min_GPa=scores["pressure_histogram_min_GPa"],
        pressure_histogram_max_GPa=scores["pressure_histogram_max_GPa"],
    )


def calc_md_pressure_metrics(df_eval: pd.DataFrame) -> dict[str, float]:
    """Aggregate source pressure percentage-error rows for one model."""
    status = df_eval.get("status", pd.Series(dtype=str)).fillna("")
    if "pressure_error_percent" in df_eval:
        error_percent = pd.to_numeric(
            df_eval["pressure_error_percent"], errors="coerce"
        )
    elif "final_mean_pressure_error_percent" in df_eval:
        error_percent = pd.to_numeric(
            df_eval["final_mean_pressure_error_percent"], errors="coerce"
        )
    else:
        raise ValueError(
            "Pressure metrics require a pressure_error_percent or "
            "final_mean_pressure_error_percent column"
        )
    error_percent = error_percent.dropna()
    return {
        "pressure_error_percent": (
            float(error_percent.mean()) if len(error_percent) else float("nan")
        ),
        "n_pressure_systems": len(error_percent),
        "n_pressure_failed": int(status.eq("failed").sum()),
        "n_pressure_skipped": int(status.str.startswith("skipped").sum()),
    }


def calc_md_combined_error(metrics: dict[str, float]) -> float | None:
    """Average RDF, VDOS, and pressure errors on a shared 0-to-1 scale."""
    rdf_error = metrics.get("rdf_error")
    if rdf_error is None or not np.isfinite(float(rdf_error)):
        return None
    rdf_error_fraction = max(0.0, float(rdf_error) / 100)
    vdos_error_value = metrics.get("vdos_error")
    if vdos_error_value is None or not np.isfinite(float(vdos_error_value)):
        return None
    vdos_error_fraction = max(0.0, float(vdos_error_value))
    pressure_error = metrics.get("pressure_error_percent")
    if pressure_error is None or not np.isfinite(float(pressure_error)):
        return None
    pressure_error_fraction = max(0.0, float(pressure_error) / 100)

    errors = np.asarray(
        [
            rdf_error_fraction,
            vdos_error_fraction,
            pressure_error_fraction,
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(errors)):
        return None
    return float(errors.mean())


def add_md_combined_error(metrics: dict[str, float]) -> dict[str, float]:
    """Return MD metrics with a combined error when every component is valid."""
    combined_error = calc_md_combined_error(metrics)
    if combined_error is None:
        return metrics
    return metrics | {"combined_error": combined_error}


def calc_md_metrics(df_eval: pd.DataFrame) -> dict[str, float]:
    """Aggregate a standalone MD prediction file or a bundled MD artifact."""
    aggregators = (
        (
            "reference",
            ("energy_error_per_atom", "force_sse"),
            calc_md_reference_metrics,
        ),
        ("rdf", ("rdf_error", "RDF_Error"), calc_md_rdf_metrics),
        (
            "vdos",
            ("vdos_error", "final_mean_vdos_error", "mean_vdos_error"),
            calc_md_vdos_metrics,
        ),
        (
            "pressure",
            ("pressure_error_percent", "final_mean_pressure_error_percent"),
            calc_md_pressure_metrics,
        ),
    )
    metrics: dict[str, float] = {}
    if "metric_kind" in df_eval:
        for metric_kind, _, aggregator in aggregators:
            component_rows = df_eval.loc[df_eval["metric_kind"].eq(metric_kind)]
            if len(component_rows):
                metrics |= aggregator(component_rows)
        if not metrics:
            raise ValueError("No recognized MD metric kinds found in prediction bundle")
    else:
        for _, indicator_columns, aggregator in aggregators:
            if any(column in df_eval for column in indicator_columns):
                metrics |= aggregator(df_eval)
        if not metrics:
            deprecated_md_columns = {
                "rdf_similarity",
                "rdf_similarity_error",
                "vdos_similarity",
                "vdos_similarity_error",
                "pressure_similarity",
                "pressure_similarity_percent",
                "pressure_error_fraction",
            }
            if any(column in df_eval for column in deprecated_md_columns):
                raise ValueError(
                    "MD metrics require RDF/VDOS/pressure error columns, not "
                    "legacy similarity columns"
                )
            # Preserve support for older reference-frame prediction tables.
            metrics = calc_md_reference_metrics(df_eval)
    return add_md_combined_error(metrics)


def write_md_prediction_bundle(
    *,
    model_key: str,
    output_dir: str | Path,
    metric_frames: dict[str, pd.DataFrame],
) -> Path:
    """Write component MD prediction rows into one reproducible artifact."""
    valid_metric_kinds = {"reference", "rdf", "vdos", "pressure"}
    unknown_metric_kinds = set(metric_frames) - valid_metric_kinds
    if unknown_metric_kinds:
        raise ValueError(f"Unknown MD metric kinds: {sorted(unknown_metric_kinds)}")
    if not metric_frames:
        raise ValueError("At least one MD metric table is required for a bundle")

    bundle = pd.concat(
        [
            component_frame.assign(metric_kind=metric_kind)
            for metric_kind, component_frame in metric_frames.items()
        ],
        ignore_index=True,
        sort=False,
    )
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    pred_path = output_path / f"md_predictions_bundle_{model_key}.csv"
    bundle.to_csv(pred_path, index=False)
    return pred_path


def evaluate_pressure_dataset(
    *,
    model_key: str,
    calculator: Calculator,
    ref_dir: str | Path,
    mlip_dir: str | Path,
    output_dir: str | Path,
    config: MdPressureConfig | None = None,
    md_config: MdRunConfig | None = None,
    extensions: Sequence[str] = (".extxyz", ".xyz", ".extxyz.gz", ".xyz.gz"),
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Evaluate matched-length pressure histogram error for MLIP rollouts."""
    config = config or MdPressureConfig()
    md_config = md_config or MdRunConfig()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    ref_settings = load_rdf_time_step_settings(config.ref_settings_path)
    mlip_settings = load_rdf_time_step_settings(config.mlip_settings_path)
    rows: list[dict[str, Any]] = []
    all_ref_files = find_trajectory_files(ref_dir, extensions=extensions)

    with tqdm(
        total=len(all_ref_files),
        desc=f"Evaluating pressure with {model_key}",
        unit="system",
    ) as progress_bar:
        for ref_file in all_ref_files:
            system_name = ref_file.parent.name
            progress_bar.set_postfix_str(system_name, refresh=False)
            skip_pressure = system_name.startswith(config.skip_dir_prefixes)
            if skip_pressure or should_skip_trajectory(ref_file, md_config):
                status = "skipped_pressure_policy" if skip_pressure else "skipped_input"
                rows.append(
                    asdict(
                        MdPressureEvalResult(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=str(ref_file),
                            mlip_file=None,
                            pressure_error_percent=None,
                            pressure_histogram_l1_area=None,
                            pressure_histogram_distance=None,
                            pressure_ref_histogram_area=None,
                            pressure_mlip_histogram_area=None,
                            pressure_histogram_min_GPa=None,
                            pressure_histogram_max_GPa=None,
                            bins=config.bins,
                            n_ref_frames=0,
                            n_mlip_frames=0,
                            ref_frame_dt_fs=None,
                            mlip_frame_dt_fs=None,
                            matched_time_fs=None,
                            status=status,
                        )
                    )
                )
                progress_bar.update()
                continue

            mlip_file = find_mlip_trajectory(
                mlip_dir=mlip_dir,
                system_name=system_name,
                model_key=model_key,
            )
            try:
                if mlip_file is None:
                    raise FileNotFoundError(f"No MLIP trajectory found in {mlip_dir}")
                ref_dt_fs = config.ref_frame_dt_fs or _frame_dt_from_settings(
                    system_name, ref_settings
                )
                mlip_dt_fs = config.mlip_frame_dt_fs or _frame_dt_from_settings(
                    system_name, mlip_settings
                )
                if ref_dt_fs is None or mlip_dt_fs is None:
                    raise ValueError(
                        "Pressure evaluation needs reference and MLIP frame "
                        "spacings or matching settings CSVs"
                    )
                rows.append(
                    asdict(
                        evaluate_pressure_trajectory_pair(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=ref_file,
                            mlip_file=mlip_file,
                            output_dir=output_path,
                            calculator=calculator,
                            config=config,
                            ref_frame_dt_fs=ref_dt_fs,
                            mlip_frame_dt_fs=mlip_dt_fs,
                        )
                    )
                )
            except Exception as exc:  # noqa: BLE001 - record per-system failure
                rows.append(
                    asdict(
                        MdPressureEvalResult(
                            model_key=model_key,
                            system_name=system_name,
                            ref_file=str(ref_file),
                            mlip_file=str(mlip_file) if mlip_file else None,
                            pressure_error_percent=None,
                            pressure_histogram_l1_area=None,
                            pressure_histogram_distance=None,
                            pressure_ref_histogram_area=None,
                            pressure_mlip_histogram_area=None,
                            pressure_histogram_min_GPa=None,
                            pressure_histogram_max_GPa=None,
                            bins=config.bins,
                            n_ref_frames=0,
                            n_mlip_frames=0,
                            ref_frame_dt_fs=None,
                            mlip_frame_dt_fs=None,
                            matched_time_fs=None,
                            status="failed",
                            error=f"{exc!r}\n{traceback.format_exc()}",
                        )
                    )
                )
            progress_bar.update()

    df_eval = pd.DataFrame(rows)
    pred_path = output_path / f"md_pressure_same_simulation_length_{model_key}.csv"
    df_eval.to_csv(pred_path, index=False)
    metrics = calc_md_pressure_metrics(df_eval)
    with open(
        output_path / f"md_pressure_metrics_same_simulation_length_{model_key}.json",
        mode="w",
        encoding="utf-8",
    ) as file:
        json.dump(metrics, file, indent=2)
    return df_eval, metrics


def evaluate_reference_dataset(
    *,
    model_key: str,
    calculator: Calculator,
    input_dir: str | Path,
    output_dir: str | Path,
    config: MdRunConfig | None = None,
    extensions: Sequence[str] = (".extxyz", ".xyz", ".extxyz.gz", ".xyz.gz"),
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Evaluate a model on all reference MD frames in an input directory."""
    config = config or MdRunConfig()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    rows: list[pd.DataFrame] = []
    all_trajectory_files = find_trajectory_files(input_dir, extensions=extensions)
    trajectory_files: list[Path] = []
    for file_path in all_trajectory_files:
        system_name = file_path.parent.name
        if should_skip_trajectory(file_path, config):
            rows.append(
                pd.DataFrame(
                    [
                        {
                            "model_key": model_key,
                            "input_file": str(file_path),
                            "system_name": system_name,
                            "status": "skipped_input",
                        }
                    ]
                )
            )
        else:
            trajectory_files.append(file_path)

    n_frames_total = sum(
        count_trajectory_frames(file_path) for file_path in trajectory_files
    )

    with tqdm(
        total=n_frames_total,
        desc=f"Evaluating MD refs with {model_key}",
        unit="frame",
    ) as progress_bar:
        for file_path in trajectory_files:
            system_name = file_path.parent.name

            try:
                frames = read_trajectory(file_path)
                rows.append(
                    evaluate_reference_frames(
                        frames,
                        calculator=calculator,
                        model_key=model_key,
                        input_file=file_path,
                        system_name=system_name,
                        progress_bar=progress_bar,
                    )
                )
            except Exception as exc:  # noqa: BLE001 - record per-trajectory failure
                rows.append(
                    pd.DataFrame(
                        [
                            {
                                "model_key": model_key,
                                "input_file": str(file_path),
                                "system_name": system_name,
                                "status": "failed",
                                "error": repr(exc),
                                "traceback": traceback.format_exc(),
                            }
                        ]
                    )
                )

    df_eval = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    pred_path = output_path / f"md_reference_{model_key}.csv"
    df_eval.to_csv(pred_path, index=False)

    metrics = calc_md_reference_metrics(df_eval)
    with open(
        output_path / f"md_reference_metrics_{model_key}.json",
        mode="w",
        encoding="utf-8",
    ) as file:
        json.dump(metrics, file, indent=2)

    return df_eval, metrics


def _repo_relative_path(file_path: str | Path) -> str:
    """Return a path relative to the repo root when possible."""
    path = Path(file_path)
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def _yaml_number(value: float) -> float:
    """Round finite floats while preserving NaN values for missing metrics."""
    if isinstance(value, int | np.integer):
        return int(value)
    value = float(value)
    if np.isfinite(value):
        return float(round(value, 4))
    return float("nan")


def write_metrics_to_yaml(
    model: Model,
    metrics: dict[str, float],
    pred_file_path: str | Path | None = None,
    pred_file_url: str | None = None,
) -> CommentedMap:
    """Write MD reference metrics to a model YAML file."""
    metrics = add_md_combined_error(metrics)
    deprecated_yaml_metric_keys = {
        "rdf_similarity",
        "rdf_similarity_error",
        "vdos_similarity",
        "vdos_similarity_error",
        "pressure_similarity",
        "pressure_similarity_percent",
        "pressure_error_fraction",
        "combined_similarity",
    }
    metrics_for_yaml = CommentedMap(
        {
            key: _yaml_number(value)
            for key, value in metrics.items()
            if isinstance(value, int | float | np.integer | np.floating)
            and key not in deprecated_yaml_metric_keys
        }
    )
    if pred_file_path is not None:
        metrics_for_yaml["pred_file"] = _repo_relative_path(pred_file_path)
    if pred_file_url is not None:
        metrics_for_yaml["pred_file_url"] = pred_file_url

    metric_units = {
        "energy_rmse": "eV/atom",
        "force_rmse": "eV/Angstrom",
        "rdf_error": "%",
        "vdos_error": "fraction",
        "pressure_error_percent": "%",
        "combined_error": "fraction",
        "n_frames": "count",
        "n_energy_frames": "count",
        "n_force_frames": "count",
        "n_rdf_systems": "count",
        "n_rdf_failed": "count",
        "n_rdf_skipped": "count",
        "n_vdos_systems": "count",
        "n_vdos_failed": "count",
        "n_vdos_skipped": "count",
        "n_pressure_systems": "count",
        "n_pressure_failed": "count",
        "n_pressure_skipped": "count",
        "n_failed": "count",
        "n_skipped": "count",
    }
    for key in metrics_for_yaml:
        if unit := metric_units.get(key):
            metrics_for_yaml.yaml_add_eol_comment(unit, key, column=1)

    update_yaml_file(model.yaml_path, "metrics.md", metrics_for_yaml)
    return metrics_for_yaml
