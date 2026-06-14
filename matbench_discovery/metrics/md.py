"""Metrics comparing MLIP molecular dynamics trajectories to ab-initio references.

Implements the finite-temperature observables proposed for the MD benchmark task
(https://github.com/janosh/matbench-discovery/pull/344): energy/force RMSE on
reference MD frames, radial distribution function (RDF) error, pressure MAE from
the stress tensor trace, vibrational density of states (VDOS) error from the
velocity autocorrelation function, and an error-weighted combined RDF+VDOS score.
Also adds a Wasserstein-1 distance between pressure distributions as a
frame-pairing-independent complement to the pressure MAE.
"""

import itertools
import math
from fractions import Fraction
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from ase import Atoms, units
from ase.calculators.calculator import PropertyNotImplementedError
from ase.geometry import find_mic
from ruamel.yaml.comments import CommentedMap

from matbench_discovery import ROOT
from matbench_discovery.data import update_yaml_file

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ase.calculators.calculator import Calculator

    from matbench_discovery.enums import Model

# units of per-system metric columns shared between eval scripts and aggregation
METRIC_UNITS = {
    "energy_rmse": "eV/atom",
    "force_rmse": "eV/Å",
    "rdf_error": "%",
    "vdos_error": "%",
    "pressure_mae": "GPa",
    "pressure_wasserstein": "GPa",
    "combined_error": "%",  # error-weighted RDF+VDOS combination
    "n_systems": "count",
}
PER_SYSTEM_METRIC_COLS = tuple(
    key for key in METRIC_UNITS if key not in ("combined_error", "n_systems")
)


def calc_rdf(
    trajectory: "Sequence[Atoms]", *, n_bins: int = 500, r_max: float | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Time-averaged all-pair radial distribution function g(r) of a trajectory.

    Returns bin-center distances r in Angstrom and g(r), histogrammed into n_bins
    between 0 and r_max. r_max defaults to half the smallest cell length of the
    first frame, the validity limit of minimum-image distances.
    """
    if len(trajectory) == 0:
        raise ValueError("Cannot compute RDF of empty trajectory")
    if n_bins < 2:
        raise ValueError(f"Need >= 2 RDF bins, got {n_bins=}")
    n_atoms = len(trajectory[0])
    if n_atoms < 2:
        raise ValueError(f"Cannot compute RDF with {n_atoms=} < 2")
    min_half_cell = min(float(min(atoms.cell.lengths())) for atoms in trajectory) / 2
    if r_max is None:
        r_max = min_half_cell
    elif r_max > min_half_cell + 1e-12:
        raise ValueError(
            f"{r_max=} exceeds minimum-image validity limit {min_half_cell:.6g} A "
            "for at least one trajectory frame"
        )
    if r_max <= 0:
        raise ValueError(f"{r_max=} must be positive")

    bin_edges = np.linspace(0, r_max, n_bins + 1)
    counts = np.zeros(n_bins)
    upper_tri_idx = np.triu_indices(n_atoms, k=1)
    inv_volume_sum = 0.0  # per-frame densities make this exact for varying cells
    for atoms in trajectory:
        if len(atoms) != n_atoms:
            raise ValueError(
                f"Inconsistent atom counts in trajectory: {len(atoms)=} != {n_atoms=}"
            )
        pair_dists = atoms.get_all_distances(mic=True)[upper_tri_idx]
        counts += np.histogram(pair_dists, bins=bin_edges)[0]
        inv_volume_sum += 1 / atoms.get_volume()

    radii = (bin_edges[:-1] + bin_edges[1:]) / 2
    shell_volumes = 4 / 3 * np.pi * np.diff(bin_edges**3)
    # Paper Eq. 4 normalization: rho*N shell volume. The upper-triangle histogram
    # counts unordered pairs, hence the /2 relative to Eq. 4's ordered i != j sum.
    ideal_counts = n_atoms**2 / 2 * inv_volume_sum * shell_volumes
    return radii, counts / ideal_counts


def calc_rdf_error(
    radii: np.ndarray, g_r_ref: np.ndarray, g_r_pred: np.ndarray
) -> float:
    """RDF error in percent (paper Eq. 3) between two g(r) on a shared radii grid.

        error = 100 * min(1, integral |g_ref - g_pred| dr / integral |g_ref - 1| dr)

    0 means perfect match, 100 means as different from the reference as an ideal
    gas (also returned if the reference itself is indistinguishable from an ideal
    gas, i.e. zero denominator).
    """
    if not len(radii) == len(g_r_ref) == len(g_r_pred):
        raise ValueError(f"{len(radii)=}, {len(g_r_ref)=}, {len(g_r_pred)=} differ")
    g_r_ref = np.asarray(g_r_ref)
    denominator = np.trapezoid(np.abs(g_r_ref - 1), x=radii)
    if denominator == 0:
        return 100.0
    numerator = np.trapezoid(np.abs(g_r_ref - g_r_pred), x=radii)
    return float(min(1, numerator / denominator) * 100)


def calc_velocities(
    trajectory: "Sequence[Atoms]", *, time_step_fs: float
) -> np.ndarray:
    """Finite-difference velocities of shape (n_frames, n_atoms, 3) in Angstrom/fs.

    Positions are unwrapped via minimum-image displacements between consecutive
    frames before differentiating, so periodic boundary crossings don't produce
    velocity spikes. Using the same estimator for reference and MLIP trajectories
    keeps VDOS comparisons unbiased even when one stores explicit velocities.
    """
    if len(trajectory) < 2:
        raise ValueError(
            f"Need >= 2 frames to estimate velocities, got {len(trajectory)}"
        )
    if time_step_fs <= 0:
        raise ValueError(f"{time_step_fs=} must be positive")

    mic_displacements = [
        find_mic(frame_next.positions - frame.positions, frame.cell, frame.pbc)[0]
        for frame, frame_next in itertools.pairwise(trajectory)
    ]
    unwrapped = trajectory[0].positions + np.cumsum(
        [np.zeros_like(trajectory[0].positions), *mic_displacements], axis=0
    )
    return np.gradient(unwrapped, time_step_fs, axis=0)


def calc_vdos(
    velocities: np.ndarray, *, time_step_fs: float
) -> tuple[np.ndarray, np.ndarray]:
    """Vibrational density of states from velocities of shape (n_frames, n_atoms, 3).

    Returns frequencies in THz and VDOS intensities: the one-sided power spectrum of
    each Hann-windowed velocity component (the Fourier transform of its
    autocorrelation function by the Wiener-Khinchin theorem), normalized to unit
    total power per component before summing. Per-component normalization weights
    every degree of freedom equally as required by equipartition (a raw sum would
    overweight light atoms by 1/mass) and mirrors the paper's normalized-VACF
    estimator. Components with zero total power (frozen atoms) are skipped.
    """
    velocities = np.asarray(velocities, dtype=float)
    if velocities.ndim != 3 or velocities.shape[-1] != 3:
        raise ValueError(
            f"Expected shape (n_frames, n_atoms, 3), got {velocities.shape}"
        )
    # hanning(n < 4) has at most one non-zero sample, making the windowed
    # spectrum identically zero or a featureless constant
    if (n_frames := len(velocities)) < 4:
        raise ValueError(f"Need >= 4 frames to compute VDOS, got {n_frames}")
    if time_step_fs <= 0:
        raise ValueError(f"{time_step_fs=} must be positive")

    window = np.hanning(n_frames)[:, None, None]
    spectra = np.abs(np.fft.rfft(velocities * window, axis=0)) ** 2
    total_power = spectra.sum(axis=0, keepdims=True)
    spectra = np.divide(
        spectra, total_power, out=np.zeros_like(spectra), where=total_power > 0
    )
    freqs_thz = np.fft.rfftfreq(n_frames, d=time_step_fs * 1e-3)  # fs -> ps = 1/THz
    return freqs_thz, spectra.sum(axis=(1, 2))


def calc_vdos_error(
    freqs_ref: np.ndarray,
    vdos_ref: np.ndarray,
    freqs_pred: np.ndarray,
    vdos_pred: np.ndarray,
) -> float:
    """VDOS error in percent (paper Eq. 7): 0 for a perfect match, 100 for
    non-overlapping spectra. Both spectra are area-normalized and interpolated onto
    the union frequency grid up to the smaller of the two Nyquist frequencies.

        error = 100 * integral |f_ref - f_pred| / (integral f_ref + integral f_pred)
    """
    f_max = min(np.max(freqs_ref), np.max(freqs_pred))
    grid = np.union1d(freqs_ref, freqs_pred)
    grid = grid[grid <= f_max]
    if len(grid) < 2:
        raise ValueError("Frequency grids have insufficient overlap")

    ref_interp = np.interp(grid, freqs_ref, vdos_ref)
    pred_interp = np.interp(grid, freqs_pred, vdos_pred)
    for label, spectrum in (("reference", ref_interp), ("predicted", pred_interp)):
        if (spectrum < 0).any():
            raise ValueError(f"{label} VDOS has negative intensities")
        if np.trapezoid(spectrum, x=grid) <= 0:
            raise ValueError(f"{label} VDOS has non-positive area")

    ref_interp = ref_interp / np.trapezoid(ref_interp, x=grid)
    pred_interp = pred_interp / np.trapezoid(pred_interp, x=grid)
    numerator = np.trapezoid(np.abs(ref_interp - pred_interp), x=grid)
    # both spectra are area-normalized, so Eq. 7's denominator (sum of areas) == 2
    return float(numerator / 2 * 100)


def calc_pressure(stress: np.ndarray) -> float:
    """Pressure in GPa from an ASE stress tensor (paper Eq. 5): P = -Tr(sigma)/3.

    Accepts stress in eV/Å³ as Voigt 6-vector [xx, yy, zz, yz, xz, xy] or (3, 3)
    matrix, following ASE conventions.
    """
    stress = np.asarray(stress, dtype=float)
    if stress.shape == (6,):
        trace = stress[:3].sum()
    elif stress.shape == (3, 3):
        trace = np.trace(stress)
    else:
        raise ValueError(f"Expected stress of shape (6,) or (3, 3), got {stress.shape}")
    return float(-trace / 3 / units.GPa)


def get_trajectory_pressures(trajectory: "Sequence[Atoms]") -> np.ndarray:
    """Per-frame pressures in GPa from stresses stored on trajectory frames (e.g.
    attached via SinglePointCalculator when reading extxyz files). Raises
    RuntimeError or PropertyNotImplementedError if any frame lacks stress.
    """
    return np.array(
        [calc_pressure(atoms.get_stress(voigt=True)) for atoms in trajectory]
    )


def _time_strides(ref_time_step_fs: float, pred_time_step_fs: float) -> tuple[int, int]:
    """Frame strides (ref, pred) that subsample two trajectories with different
    saved-frame intervals onto their common (least common multiple) time grid, i.e.
    frames i * ref_stride and i * pred_stride share identical timestamps.
    """
    if min(ref_time_step_fs, pred_time_step_fs) <= 0:
        raise ValueError(
            f"Time steps must be positive, got {ref_time_step_fs=}, "
            f"{pred_time_step_fs=}"
        )
    ref_time_step = Fraction(str(ref_time_step_fs))
    pred_time_step = Fraction(str(pred_time_step_fs))
    common_time_step = Fraction(
        math.lcm(ref_time_step.numerator, pred_time_step.numerator),
        math.gcd(ref_time_step.denominator, pred_time_step.denominator),
    )
    return int(common_time_step / ref_time_step), int(common_time_step / pred_time_step)


def calc_pressure_metrics(
    pressures_ref: np.ndarray,
    pressures_pred: np.ndarray,
    *,
    ref_time_step_fs: float = 1,
    pred_time_step_fs: float = 1,
) -> dict[str, float]:
    """Pressure MAE (paper Eq. 6) and Wasserstein-1 distance in GPa.

    The MAE pairs frames at identical timestamps: both arrays are subsampled onto
    their common (LCM) time grid before pairing, so different saved-frame cadences
    compare the same simulation times. It still conflates mean pressure offsets
    with fluctuation decorrelation between independently thermalized runs. The W1
    distance compares the full pressure distributions instead, making it
    insensitive to frame pairing while still capturing mean offsets and
    fluctuation-width differences.
    """
    from scipy.stats import wasserstein_distance

    pressures_ref, pressures_pred = map(np.asarray, (pressures_ref, pressures_pred))
    if len(pressures_ref) == 0 or len(pressures_pred) == 0:
        raise ValueError("Cannot compute pressure metrics of empty pressure arrays")

    ref_stride, pred_stride = _time_strides(ref_time_step_fs, pred_time_step_fs)
    ref_aligned = pressures_ref[::ref_stride]
    pred_aligned = pressures_pred[::pred_stride]
    n_pairs = min(len(ref_aligned), len(pred_aligned))
    mae = np.abs(ref_aligned[:n_pairs] - pred_aligned[:n_pairs]).mean()
    return {
        "pressure_mae": float(mae),
        "pressure_wasserstein": float(
            wasserstein_distance(pressures_ref, pressures_pred)
        ),
    }


def calc_combined_error(rdf_error: float, vdos_error: float) -> float:
    """Error complement of the paper Eq. 8 combined RDF/VDOS score in percent.

        combined_score = eps * (100 - rdf_error) + (1 - eps) * (100 - vdos_error)
        eps = rdf_error / (rdf_error + vdos_error)
        combined_error = 100 - combined_score

    This preserves the leaderboard's lower-is-better error convention while using
    Eq. 8's weighting: poor performance receives more weight, so a model can't hide
    a catastrophic VDOS behind a good RDF.
    """
    if rdf_error < 0 or vdos_error < 0:
        raise ValueError(
            f"Errors must be non-negative, got {rdf_error=}, {vdos_error=}"
        )
    if (total := rdf_error + vdos_error) == 0:
        return 0.0
    # eps * rdf + (1 - eps) * vdos simplifies to (rdf^2 + vdos^2) / (rdf + vdos)
    return float((rdf_error**2 + vdos_error**2) / total)


def calc_energy_force_rmse(
    trajectory: "Sequence[Atoms]", calculator: "Calculator"
) -> dict[str, float]:
    """Energy (eV/atom) and force (eV/Å) RMSE of calculator single-point predictions
    on reference MD frames carrying DFT energies and forces (paper Eq. 1-2).

    Energies are compared per atom, preserving constant offsets as specified by Eq. 1.
    This means energy RMSE can be dominated by differing absolute energy references.
    """
    if len(trajectory) == 0:
        raise ValueError("Cannot compute energy/force RMSE of empty trajectory")

    energy_errs_per_atom: list[float] = []
    force_sq_err_sum, n_force_components = 0.0, 0
    for atoms in trajectory:
        ref_energy, ref_forces = atoms.get_potential_energy(), atoms.get_forces()
        pred_atoms = atoms.copy()
        pred_atoms.calc = calculator
        energy_errs_per_atom.append(
            (pred_atoms.get_potential_energy() - ref_energy) / len(atoms)
        )
        force_sq_err_sum += ((pred_atoms.get_forces() - ref_forces) ** 2).sum()
        n_force_components += ref_forces.size

    return {
        "energy_rmse": float(np.sqrt(np.mean(np.square(energy_errs_per_atom)))),
        "force_rmse": float(np.sqrt(force_sq_err_sum / n_force_components)),
    }


def matched_frame_counts(
    *,
    n_ref_frames: int,
    n_pred_frames: int,
    ref_time_step_fs: float,
    pred_time_step_fs: float,
) -> tuple[int, int]:
    """Frame counts (n_ref_use, n_pred_use) covering the identical closed time
    interval [0, t_match] in both trajectories, where t_match is the largest
    multiple of the common (LCM) time grid that fits both trajectory durations.
    A trajectory of n frames saved every dt spans (n - 1) * dt.
    """
    ref_stride, pred_stride = _time_strides(ref_time_step_fs, pred_time_step_fs)
    if min(n_ref_frames, n_pred_frames) < 1:
        return 0, 0
    # number of common-grid steps fitting both trajectory durations
    n_common = min((n_ref_frames - 1) // ref_stride, (n_pred_frames - 1) // pred_stride)
    # + 1 to include the frame at t=0
    return n_common * ref_stride + 1, n_common * pred_stride + 1


def evaluate_md_system(
    ref_trajectory: "Sequence[Atoms]",
    pred_trajectory: "Sequence[Atoms]",
    *,
    ref_time_step_fs: float,
    pred_time_step_fs: float,
    calculator: "Calculator | None" = None,
    n_rdf_bins: int = 500,
) -> dict[str, float]:
    """All per-system MD metrics for one reference/predicted trajectory pair, keyed
    by METRIC_UNITS column names.

    Both trajectories are truncated to the same total simulation time (given their
    respective time steps in fs between saved frames) before computing observables,
    matching the paper's same-simulation-length protocol. Pressure metrics are NaN
    when either trajectory lacks stress data. Energy and force RMSE on the full
    reference trajectory are included when a calculator is passed.
    """
    n_ref_use, n_pred_use = matched_frame_counts(
        n_ref_frames=len(ref_trajectory),
        n_pred_frames=len(pred_trajectory),
        ref_time_step_fs=ref_time_step_fs,
        pred_time_step_fs=pred_time_step_fs,
    )
    if n_ref_use < 4 or n_pred_use < 4:  # VDOS needs >= 4 frames
        raise ValueError(
            f"Trajectories too short after time matching: {n_ref_use=}, {n_pred_use=}"
        )
    ref_matched = ref_trajectory[:n_ref_use]
    pred_matched = pred_trajectory[:n_pred_use]

    radii, g_r_ref = calc_rdf(ref_matched, n_bins=n_rdf_bins)
    r_max = float(radii[-1] + (radii[1] - radii[0]) / 2)
    _, g_r_pred = calc_rdf(pred_matched, n_bins=n_rdf_bins, r_max=r_max)
    metrics = {"rdf_error": calc_rdf_error(radii, g_r_ref, g_r_pred)}

    ref_vels = calc_velocities(ref_matched, time_step_fs=ref_time_step_fs)
    pred_vels = calc_velocities(pred_matched, time_step_fs=pred_time_step_fs)
    freqs_ref, vdos_ref = calc_vdos(ref_vels, time_step_fs=ref_time_step_fs)
    freqs_pred, vdos_pred = calc_vdos(pred_vels, time_step_fs=pred_time_step_fs)
    metrics["vdos_error"] = calc_vdos_error(freqs_ref, vdos_ref, freqs_pred, vdos_pred)

    try:
        metrics |= calc_pressure_metrics(
            get_trajectory_pressures(ref_matched),
            get_trajectory_pressures(pred_matched),
            ref_time_step_fs=ref_time_step_fs,
            pred_time_step_fs=pred_time_step_fs,
        )
    except (PropertyNotImplementedError, RuntimeError) as exc:
        # NaN pressure metrics only for frames without (a calculator providing)
        # stress. NB: PropertyNotImplementedError subclasses RuntimeError via
        # NotImplementedError, so check it first; other RuntimeErrors re-raise.
        if not isinstance(exc, PropertyNotImplementedError) and (
            "no calculator" not in str(exc)
        ):
            raise
        metrics["pressure_mae"] = metrics["pressure_wasserstein"] = float("nan")

    if calculator is not None:
        metrics |= calc_energy_force_rmse(ref_trajectory, calculator)

    return metrics


def load_per_system_metrics(csv_paths: "Sequence[str]") -> pd.DataFrame:
    """Concatenate per-system MD metric CSVs (one row each, as written by parallel
    single-system runs) into a single DataFrame indexed by system.

    Deduplicates on the system column keeping the last occurrence, so re-running a
    system (e.g. a timed-out rollout finished in a later job) overrides the earlier
    row instead of double-counting it.
    """
    if not csv_paths:
        raise ValueError("No per-system metric CSVs given")
    df_all = pd.concat([pd.read_csv(path) for path in csv_paths], ignore_index=True)
    if "system" not in df_all:
        raise ValueError(f"per-system CSVs lack a 'system' column, got {[*df_all]}")
    return df_all.drop_duplicates(subset="system", keep="last").set_index("system")


def calc_md_metrics(df_md: pd.DataFrame) -> dict[str, float]:
    """Aggregate per-system MD metric rows (one row per system with a subset of
    PER_SYSTEM_METRIC_COLS columns) into model-level metrics.

    Means are taken over systems, skipping NaNs (e.g. systems without stress data).
    The combined RDF+VDOS error is the lower-is-better complement of the paper's
    Eq. 8 Pareto score, computed from model-level mean errors.
    """
    metric_cols = [col for col in PER_SYSTEM_METRIC_COLS if col in df_md]
    if not metric_cols:
        raise ValueError(
            f"No recognized MD metric columns in {list(df_md.columns)=}, "
            f"expected a subset of {PER_SYSTEM_METRIC_COLS}"
        )
    metrics = {col: float(df_md[col].mean()) for col in metric_cols}
    if "rdf_error" in metrics and "vdos_error" in metrics:
        metrics["combined_error"] = calc_combined_error(
            metrics["rdf_error"], metrics["vdos_error"]
        )
    metrics["n_systems"] = len(df_md)
    return metrics


def write_metrics_to_yaml(
    model: "Model",
    metrics: dict[str, float],
    pred_file_path: str | None = None,
    pred_file_url: str | None = None,
) -> dict[str, Any]:
    """Write model-level MD metrics from calc_md_metrics to the metrics.md section
    of a model YAML file, with unit comments and optional pred_file references.
    """
    yaml_metrics = CommentedMap()
    if pred_file_path is not None:
        yaml_metrics["pred_file"] = pred_file_path.removeprefix(f"{ROOT}/")
    if pred_file_url is not None:
        yaml_metrics["pred_file_url"] = pred_file_url
    for key, value in metrics.items():
        yaml_metrics[key] = int(value) if key == "n_systems" else float(round(value, 4))
        if unit := METRIC_UNITS.get(key):
            yaml_metrics.yaml_add_eol_comment(unit, key, column=1)

    update_yaml_file(model.yaml_path, "metrics.md", yaml_metrics)
    return yaml_metrics
