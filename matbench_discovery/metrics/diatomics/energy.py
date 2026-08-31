"""Energy-based metrics for diatomic curves."""

import numpy as np
from ase.data import atomic_masses, atomic_numbers
from numpy.typing import ArrayLike

PBE_WALL_ENERGY_THRESHOLDS_EV: tuple[float, ...] = (1, 5, 10, 20, 50, 100)


def _validate_diatomic_curve(
    xs: ArrayLike, ys: ArrayLike
) -> tuple[np.ndarray, np.ndarray]:
    """Validate curve input data.

    Args:
        xs (ArrayLike[float]): interatomic distances
        ys (ArrayLike[Any]): Energies or forces

    Returns:
        tuple[np.ndarray, np.ndarray]: Validated x and y arrays sorted by ascending x

    Raises:
        ValueError: If input data is invalid
    """
    xs_arr: np.ndarray = np.asarray(xs)
    ys_arr: np.ndarray = np.asarray(ys)

    if len(xs_arr) != len(ys_arr):
        raise ValueError(f"{len(xs_arr)=} != {len(ys_arr)=}")
    if len(xs_arr) < 2:
        raise ValueError(f"Input must have at least 2 points, got {len(xs_arr)=}")
    n_x_nan, n_y_nan = int(np.isnan(xs_arr).sum()), int(np.isnan(ys_arr).sum())
    if n_x_nan or n_y_nan:
        raise ValueError(f"Input contains NaN values: {n_x_nan=}, {n_y_nan=}")
    n_x_inf, n_y_inf = int(np.isinf(xs_arr).sum()), int(np.isinf(ys_arr).sum())
    if n_x_inf or n_y_inf:
        raise ValueError(f"Input contains infinite values: {n_x_inf=}, {n_y_inf=}")
    n_unique = len(np.unique(xs_arr))
    if n_unique != len(xs_arr):
        raise ValueError(f"xs contains {len(xs_arr) - n_unique} duplicates")

    sort_idx = np.argsort(xs_arr)  # ascending order
    return xs_arr[sort_idx], ys_arr[sort_idx]


def _validate_curve_pair(
    seps_ref: ArrayLike,
    ys_ref: ArrayLike,
    seps_pred: ArrayLike,
    ys_pred: ArrayLike,
    *,
    interpolate: bool | int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Validate/sort a ref/pred curve pair, requiring same x-grid if not
    interpolating. Returns sorted (seps_ref, ys_ref, seps_pred, ys_pred) arrays.
    """
    seps_ref, ys_ref = _validate_diatomic_curve(seps_ref, ys_ref)
    seps_pred, ys_pred = _validate_diatomic_curve(seps_pred, ys_pred)
    if not interpolate and not np.array_equal(seps_ref, seps_pred):
        raise ValueError(
            f"Reference and predicted distances must be same when {interpolate=}\n"
            f"{seps_ref=}, {seps_pred=}"
        )
    return seps_ref, ys_ref, seps_pred, ys_pred


def _interp_curve(seps_new: np.ndarray, seps: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Linearly interpolate a curve onto new distances, one trailing component at a
    time so energies (n_seps,) and forces (n_seps, n_atoms, 3) share one code path.
    """
    ys_columns = ys.reshape(len(seps), -1)
    ys_interp = np.stack(
        [np.interp(seps_new, seps, column) for column in ys_columns.T], axis=1
    )
    return ys_interp.reshape(len(seps_new), *ys.shape[1:])


def _common_grid_pair(
    seps_ref: ArrayLike,
    ys_ref: ArrayLike,
    seps_pred: ArrayLike,
    ys_pred: ArrayLike,
    *,
    interpolate: bool | int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ref/pred energies or forces on their shared distance range.

    Without interpolation both curves must share one distance grid, which is returned
    as is. Otherwise both are interpolated onto a linspace of 100 (interpolate=True)
    or `interpolate` points spanning the overlap of the sampled ranges.
    """
    seps_ref, ys_ref, seps_pred, ys_pred = _validate_curve_pair(
        seps_ref, ys_ref, seps_pred, ys_pred, interpolate=interpolate
    )
    if not interpolate:
        return seps_ref, ys_ref, ys_pred

    data_min = max(seps_ref.min(), seps_pred.min())
    data_max = min(seps_ref.max(), seps_pred.max())
    # >= (not >) to also reject a single shared point: interpolating one point across
    # a grid is meaningless.
    if data_min >= data_max:
        raise ValueError(
            f"Cannot interpolate curves with no overlap: {data_min=}, {data_max=}"
        )
    n_points = 100 if interpolate is True else int(interpolate)
    seps_interp = np.linspace(data_min, data_max, n_points)
    ys_ref_interp = _interp_curve(seps_interp, seps_ref, ys_ref)
    ys_pred_interp = _interp_curve(seps_interp, seps_pred, ys_pred)
    return seps_interp, ys_ref_interp, ys_pred_interp


def _binding_energy(energies: np.ndarray) -> float:
    """Return well depth relative to the largest sampled separation."""
    return float(energies[-1] - np.min(energies))


def _quadratic_well_fit(
    seps: ArrayLike, energies: ArrayLike, n_fit_points: int = 5
) -> tuple[float, float]:
    """Estimate equilibrium separation and curvature from a local quadratic fit."""
    seps, energies = _validate_diatomic_curve(seps, energies)
    min_idx = int(np.argmin(energies))
    if len(seps) < 3:
        return float(seps[min_idx]), np.nan

    start_idx = min(
        max(0, min_idx - n_fit_points // 2), max(0, len(seps) - n_fit_points)
    )
    fit_seps = seps[start_idx : start_idx + n_fit_points]
    fit_energies = energies[start_idx : start_idx + n_fit_points]
    if len(fit_seps) < 3:
        return float(seps[min_idx]), np.nan

    quadratic_coef, linear_coef, _constant_coef = np.polyfit(fit_seps, fit_energies, 2)
    curvature = 2 * quadratic_coef
    if quadratic_coef <= 0:
        return float(seps[min_idx]), np.nan

    equilibrium_dist = -linear_coef / (2 * quadratic_coef)
    if fit_seps.min() <= equilibrium_dist <= fit_seps.max():
        return float(equilibrium_dist), float(curvature)
    return float(seps[min_idx]), float(curvature)


def _repulsive_radius_at_threshold(
    seps: ArrayLike, energies: ArrayLike, threshold_ev: float
) -> float:
    """Invert the repulsive branch to r at E_min + threshold_ev."""
    seps, energies = _validate_diatomic_curve(seps, energies)
    min_idx = int(np.argmin(energies))
    if min_idx == 0:
        return np.nan

    radii_inward = seps[min_idx::-1]
    energy_above_min = energies[min_idx::-1] - energies[min_idx]
    monotonic_energy = np.maximum.accumulate(energy_above_min)
    unique_energy, unique_idx = np.unique(monotonic_energy, return_index=True)
    if len(unique_energy) < 2 or threshold_ev > unique_energy[-1]:
        return np.nan
    return float(np.interp(threshold_ev, unique_energy, radii_inward[unique_idx]))


def calc_pbe_wall_dist_mae(
    seps_ref: ArrayLike,
    energy_ref: ArrayLike,
    seps_pred: ArrayLike,
    energy_pred: ArrayLike,
    *,
    thresholds_ev: tuple[float, ...] = PBE_WALL_ENERGY_THRESHOLDS_EV,
) -> float:
    """Mean wall-radius error from 1 to 100 eV above the well.

    Thresholds beyond the reference wall are skipped. If the reference reaches a
    threshold but the prediction does not within the scored range, the missing crossing
    receives the full reference-radius error instead of being silently dropped.
    """
    errors = []
    for threshold_ev in thresholds_ev:
        radius_ref = _repulsive_radius_at_threshold(seps_ref, energy_ref, threshold_ev)
        if not np.isfinite(radius_ref):
            continue
        radius_pred = _repulsive_radius_at_threshold(
            seps_pred, energy_pred, threshold_ev
        )
        errors.append(
            abs(radius_pred - radius_ref) if np.isfinite(radius_pred) else radius_ref
        )
    return float(np.mean(errors)) if errors else np.nan


def calc_pbe_energy_mae(
    seps_ref: ArrayLike,
    energy_ref: ArrayLike,
    seps_pred: ArrayLike,
    energy_pred: ArrayLike,
    *,
    interpolate: bool | int = 200,
) -> float:
    """Mean absolute PBE energy error after aligning both curves at dissociation."""
    _, energy_ref, energy_pred = _common_grid_pair(
        seps_ref, energy_ref, seps_pred, energy_pred, interpolate=interpolate
    )
    energy_ref = energy_ref - energy_ref[-1]
    energy_pred = energy_pred - energy_pred[-1]
    return float(np.mean(np.abs(energy_pred - energy_ref)))


def calc_pbe_bond_length_error(
    seps_ref: ArrayLike,
    energy_ref: ArrayLike,
    seps_pred: ArrayLike,
    energy_pred: ArrayLike,
    *,
    min_ref_binding_ev: float = 0.05,
) -> float:
    """Absolute equilibrium bond-length error relative to PBE."""
    seps_ref, energy_ref = _validate_diatomic_curve(seps_ref, energy_ref)
    seps_pred, energy_pred = _validate_diatomic_curve(seps_pred, energy_pred)
    if _binding_energy(energy_ref) < min_ref_binding_ev:
        return np.nan
    ref_dist = _quadratic_well_fit(seps_ref, energy_ref)[0]
    pred_dist = _quadratic_well_fit(seps_pred, energy_pred)[0]
    return float(abs(pred_dist - ref_dist))


def calc_pbe_well_depth_error(
    seps_ref: ArrayLike,
    energy_ref: ArrayLike,
    seps_pred: ArrayLike,
    energy_pred: ArrayLike,
    *,
    min_ref_binding_ev: float = 0.05,
) -> float:
    """Absolute well-depth error, D_e = E(r_max) - E_min, relative to PBE."""
    _, energy_ref = _validate_diatomic_curve(seps_ref, energy_ref)
    _, energy_pred = _validate_diatomic_curve(seps_pred, energy_pred)
    ref_depth = _binding_energy(energy_ref)
    if ref_depth < min_ref_binding_ev:
        return np.nan
    return float(abs(_binding_energy(energy_pred) - ref_depth))


def _vibrational_wavenumber_cm(elem_symbol: str, curvature_ev_per_a2: float) -> float:
    """Convert homonuclear force constant to harmonic wavenumber in cm^-1."""
    if not np.isfinite(curvature_ev_per_a2) or curvature_ev_per_a2 <= 0:
        return np.nan
    atomic_symbol = elem_symbol.split("-", maxsplit=1)[0]
    reduced_mass_kg = (
        atomic_masses[atomic_numbers[atomic_symbol]] * 1.66053906660e-27 / 2
    )
    force_constant_n_per_m = curvature_ev_per_a2 * 16.02176634
    angular_freq_per_s = np.sqrt(force_constant_n_per_m / reduced_mass_kg)
    return float(angular_freq_per_s / (2 * np.pi * 2.99792458e10))


def calc_pbe_vib_freq_error(
    elem_symbol: str,
    seps_ref: ArrayLike,
    energy_ref: ArrayLike,
    seps_pred: ArrayLike,
    energy_pred: ArrayLike,
    *,
    min_ref_binding_ev: float = 0.05,
) -> float:
    """Absolute harmonic vibrational-frequency error near the PBE well."""
    seps_ref, energy_ref = _validate_diatomic_curve(seps_ref, energy_ref)
    seps_pred, energy_pred = _validate_diatomic_curve(seps_pred, energy_pred)
    if _binding_energy(energy_ref) < min_ref_binding_ev:
        return np.nan
    ref_curvature = _quadratic_well_fit(seps_ref, energy_ref)[1]
    pred_curvature = _quadratic_well_fit(seps_pred, energy_pred)[1]
    ref_wavenumber = _vibrational_wavenumber_cm(elem_symbol, ref_curvature)
    pred_wavenumber = _vibrational_wavenumber_cm(elem_symbol, pred_curvature)
    return float(abs(pred_wavenumber - ref_wavenumber))


def calc_tortuosity(seps: ArrayLike, energies: ArrayLike) -> float:
    """Calculate tortuosity of a potential energy curve as the ratio between total
    variation in energy and the sum of absolute energy differences between shortest
    separation distance r_min, equilibrium distance r_eq, and longest separation
    distance r_max. This is essentially the arc-chord ratio projected in the energy
    dimension.

    A perfect Lennard-Jones potential or any potential with a single repulsion-
    attraction transition or pure repulsion will have tortuosity equal to 1. True PECs
    may have intermediate range energy barriers, so the elemental average should be
    slightly above 1.

    Args:
        seps (ArrayLike[float]): Interatomic distances
        energies (ArrayLike[float]): Energy values

    Returns:
        float: tortuosity value (ratio of total variation to direct energy difference).
    """
    _, energies = _validate_diatomic_curve(seps, energies)

    tv_energy = np.sum(np.abs(np.diff(energies)))
    e_min = np.min(energies)
    direct_energy_diff = abs(energies[0] - e_min) + abs(energies[-1] - e_min)

    if direct_energy_diff == 0:
        return np.nan
    return float(tv_energy / direct_energy_diff)


def _threshold_diff_signs(
    vals: np.ndarray, threshold: float = 1e-3
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute thresholded diffs, their signs, and flip mask for a 1D array.

    Args:
        vals (np.ndarray): 1D array of values (energies or forces).
        threshold (float): Diffs below this magnitude are zeroed. Defaults to 1e-3.

    Returns:
        tuple: (thresholded diffs with zeros removed, their signs, boolean flip mask)
    """
    diffs = np.diff(vals)
    diffs[np.abs(diffs) < threshold] = 0
    signs = np.sign(diffs)
    mask = signs != 0
    diffs, signs = diffs[mask], signs[mask]
    flips = np.diff(signs) != 0
    return diffs, signs, flips


def calc_energy_diff_flips(seps: ArrayLike, energies: ArrayLike) -> float:
    """Calculate number of energy difference sign flips.

    Args:
        seps (ArrayLike[float]): Interatomic distances in Å.
        energies (ArrayLike[float]): Energies in eV.

    Returns:
        float: Number of energy difference sign flips.
    """
    _, energies = _validate_diatomic_curve(seps, energies)
    _, _, flips = _threshold_diff_signs(energies)
    return float(np.sum(flips))


def calc_energy_jump(seps: ArrayLike, energies: ArrayLike) -> float:
    """Calculate energy jump metric as sum of absolute energy differences at flip
    points.

    Args:
        seps (ArrayLike[float]): Interatomic distances in Å.
        energies (ArrayLike[float]): Energies in eV.

    Returns:
        float: Sum of absolute energy differences at flip points.
    """
    _, energies = _validate_diatomic_curve(seps, energies)
    diffs, _, flips = _threshold_diff_signs(energies)
    return float(np.abs(diffs[:-1][flips]).sum() + np.abs(diffs[1:][flips]).sum())
