"""Energy-based metrics for diatomic curves."""

from collections.abc import Sequence
from typing import Any

import numpy as np


def _validate_diatomic_curve(
    xs: Sequence[float], ys: Sequence[Any], *, normalize_energy: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    """Validate curve input data.

    Args:
        xs (Sequence[float]): interatomic distances
        ys (Sequence[Any]): Energies or forces
        normalize_energy (bool): Whether to shift energies to zero at largest separation
            distance (far field). Only applies when ys are energies, not forces.

    Returns:
        tuple[np.ndarray, np.ndarray]: Validated and sorted x and y arrays

    Raises:
        ValueError: If input data is invalid
    """
    xs, ys = map(np.asarray, (xs, ys))

    if len(xs) != len(ys):
        raise ValueError(f"{len(xs)=} != {len(ys)=}")
    if len(xs) < 2:
        raise ValueError(f"Input must have at least 2 points, got {len(xs)=}")
    n_x_nan, n_y_nan = int(np.isnan(xs).sum()), int(np.isnan(ys).sum())
    if n_x_nan or n_y_nan:
        raise ValueError(f"Input contains NaN values: {n_x_nan=}, {n_y_nan=}")
    n_x_inf, n_y_inf = int(np.isinf(xs).sum()), int(np.isinf(ys).sum())
    if n_x_inf or n_y_inf:
        raise ValueError(f"Input contains infinite values: {n_x_inf=}, {n_y_inf=}")
    if len(np.unique(xs)) != len(xs):
        n_x_dup = int((np.diff(xs) == 0).sum())
        raise ValueError(f"xs contains {n_x_dup} duplicates")

    sort_idx = np.argsort(xs)[::-1]
    xs = xs[sort_idx]
    ys = ys[sort_idx]

    # If these are energies (rank 1 array), normalize to zero at far field
    if normalize_energy and ys.ndim == 1:
        ys = ys - ys[0]  # shift to zero at largest separation

    return xs, ys


def calc_curve_diff_auc(
    seps_ref: Sequence[float],
    e_ref: Sequence[float],
    seps_pred: Sequence[float],
    e_pred: Sequence[float],
    *,
    seps_range: tuple[float | None, float | None] = (None, None),
    normalize: bool = True,
) -> float:
    """Calculate the absolute area under the curve of the difference between two curves.
    Handles different x-samplings by interpolating to a common grid.

    Args:
        seps_ref (Sequence[float]): Reference interatomic distances (Å)
        e_ref (Sequence[float]): Reference potential energies (eV)
        seps_pred (Sequence[float]): Predicted interatomic distances (Å)
        e_pred (Sequence[float]): Predicted potential energies (eV)
        seps_range (tuple[float | None, float | None] | None): Optional range of
            interatomic distances to consider. Can be None to auto-set based on data
            range. If tuple is None, uses intersection of both curves' x-ranges.
        normalize (bool): Whether to normalize by reference curve's bounding box area.

    Returns:
        float: Absolute area under the curve of the difference between the curves.
            If normalize=True, returns unitless value, otherwise in eV·Å.
    """
    # Validate and sort both curves
    seps_ref, e_ref = _validate_diatomic_curve(seps_ref, e_ref, normalize_energy=False)
    seps_pred, e_pred = _validate_diatomic_curve(
        seps_pred, e_pred, normalize_energy=False
    )

    # Get data range bounds
    data_min = max(seps_ref.min(), seps_pred.min())
    data_max = min(seps_ref.max(), seps_pred.max())

    seps_min, seps_max = seps_range
    # Replace None values with data bounds
    seps_min = data_min if seps_min is None else seps_min
    seps_max = data_max if seps_max is None else seps_max

    if seps_min >= seps_max:
        raise ValueError(f"Invalid range: {seps_min=} >= {seps_max=}")

    # Create a fine grid for interpolation
    seps_interp = np.linspace(seps_min, seps_max, 1000)

    # Interpolate both curves to the common grid
    e_ref_interp = np.interp(seps_interp, seps_ref, e_ref)
    e_pred_interp = np.interp(seps_interp, seps_pred, e_pred)

    # Calculate absolute difference and integrate
    diff = np.abs(e_ref_interp - e_pred_interp)
    auc = np.trapezoid(diff, seps_interp)

    if normalize:
        # Get bounding box area of reference curve
        seps_range, e_range = np.ptp(seps_ref), np.ptp(e_ref)
        box_area = seps_range * e_range
        if box_area > 0:  # If reference curve is flat, don't normalize
            auc = auc / box_area

    return float(auc)


def calc_energy_mae(
    seps_ref: Sequence[float],
    e_ref: Sequence[float],
    seps_pred: Sequence[float],
    e_pred: Sequence[float],
) -> float:
    """Calculate mean absolute error between two energy curves.
    Handles different x-samplings by interpolating to a common grid.

    Args:
        seps_ref (Sequence[float]): Reference interatomic distances (Å)
        e_ref (Sequence[float]): Reference potential energies (eV)
        seps_pred (Sequence[float]): Predicted interatomic distances (Å)
        e_pred (Sequence[float]): Predicted potential energies (eV)

    Returns:
        float: Mean absolute error between the curves (eV).
    """
    # Validate and sort both curves
    seps_ref, e_ref = _validate_diatomic_curve(seps_ref, e_ref, normalize_energy=False)
    seps_pred, e_pred = _validate_diatomic_curve(
        seps_pred, e_pred, normalize_energy=False
    )

    # Get data range bounds
    data_min = max(seps_ref.min(), seps_pred.min())
    data_max = min(seps_ref.max(), seps_pred.max())

    # Create a fine grid for interpolation
    seps_interp = np.linspace(data_min, data_max, 1000)

    # Interpolate both curves to the common grid
    e_ref_interp = np.interp(seps_interp, seps_ref, e_ref)
    e_pred_interp = np.interp(seps_interp, seps_pred, e_pred)

    # Calculate MAE
    return float(np.mean(np.abs(e_ref_interp - e_pred_interp)))


def calc_second_deriv_smoothness(
    seps: Sequence[float], energies: Sequence[float]
) -> float:
    """Calculate smoothness using RMS of second derivative (lower is smoother)."""
    seps, energies = map(np.asarray, (seps, energies))
    sort_idx = np.argsort(seps)[::-1]  # sort in descending order
    seps = seps[sort_idx]
    energies = energies[sort_idx]
    d2y = np.gradient(np.gradient(energies, seps), seps)
    return float(np.sqrt(np.mean(d2y**2)))


def calc_total_variation_smoothness(
    seps: Sequence[float], energies: Sequence[float]
) -> float:
    """Calculate smoothness using mean absolute gradient (lower is smoother)."""
    seps, energies = map(np.asarray, (seps, energies))
    sort_idx = np.argsort(seps)[::-1]  # sort in descending order
    seps = seps[sort_idx]
    energies = energies[sort_idx]
    dy = np.gradient(energies, seps)
    return float(np.log10(np.mean(np.abs(dy))))


def calc_curvature_smoothness(
    seps: Sequence[float], energies: Sequence[float]
) -> float:
    """Calculate smoothness using mean absolute curvature (lower is smoother)."""
    seps, energies = map(np.asarray, (seps, energies))
    sort_idx = np.argsort(seps)[::-1]  # sort in descending order
    seps = seps[sort_idx]
    energies = energies[sort_idx]
    dy = np.gradient(energies, seps)
    d2y = np.gradient(dy, seps)
    curvature = np.abs(d2y) / (1 + dy**2) ** 1.5
    return float(np.log10(np.mean(curvature)))


def calc_tortuosity(seps: Sequence[float], energies: Sequence[float]) -> float:
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
        seps (Sequence[float]): Interatomic distances
        energies (Sequence[float]): Energy values

    Returns:
        float: tortuosity value (ratio of total variation to direct energy difference).
    """
    # Validate and sort with energy normalization
    _, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)

    # Total variation in energy (sum of absolute differences)
    tv_energy = np.sum(np.abs(np.diff(energies)))

    # Get minimum energy and endpoint energies
    e_min = np.min(energies)  # minimum energy (equilibrium point)
    # energy at largest distance (should be 0 after normalization)
    e_first = energies[0]
    e_last = energies[-1]  # energy at shortest distance

    # Sum of energy differences from minimum to endpoints
    direct_energy_diff = abs(e_first - e_min) + abs(e_last - e_min)

    return float(tv_energy / direct_energy_diff)


def calc_energy_diff_flips(seps: Sequence[float], energies: Sequence[float]) -> float:
    """Calculate number of energy difference sign flips.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        energies (Sequence[float]): Energies in eV.

    Returns:
        float: Number of energy difference sign flips.
    """
    seps, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)

    ediff = np.diff(energies)
    ediff[np.abs(ediff) < 1e-3] = 0  # 1meV threshold
    ediff_sign = np.sign(ediff)
    mask = ediff_sign != 0
    ediff_sign = ediff_sign[mask]
    return float(np.sum(np.diff(ediff_sign) != 0))


def calc_energy_grad_norm_max(
    seps: Sequence[float], energies: Sequence[float]
) -> float:
    """Calculate maximum absolute value of energy gradient.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        energies (Sequence[float]): Energies in eV.

    Returns:
        float: Maximum absolute value of energy gradient.
    """
    seps, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)
    return float(np.max(np.abs(np.gradient(energies, seps))))


def calc_energy_jump(seps: Sequence[float], energies: Sequence[float]) -> float:
    """Calculate energy jump metric as sum of absolute energy differences at flip
    points.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        energies (Sequence[float]): Energies in eV.

    Returns:
        float: Sum of absolute energy differences at flip points.
    """
    seps, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)

    ediff = np.diff(energies)
    ediff[np.abs(ediff) < 1e-3] = 0  # 1meV threshold
    ediff_sign = np.sign(ediff)
    mask = ediff_sign != 0
    ediff = ediff[mask]
    ediff_sign = ediff_sign[mask]
    ediff_flip = np.diff(ediff_sign) != 0

    e_jump = np.abs(ediff[:-1][ediff_flip]).sum() + np.abs(ediff[1:][ediff_flip]).sum()

    return float(e_jump)
