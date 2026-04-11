"""Energy-based metrics for diatomic curves."""

import numpy as np
from numpy.typing import ArrayLike


def _validate_diatomic_curve(
    xs: ArrayLike,
    ys: ArrayLike,
    *,
    normalize_energy: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate curve input data.

    Args:
        xs (ArrayLike[float]): interatomic distances
        ys (ArrayLike[Any]): Energies or forces
        normalize_energy (bool): Whether to shift energies to zero at largest separation
            distance (far field). Only applies when ys are energies, not forces.

    Returns:
        tuple[np.ndarray, np.ndarray]: Validated and sorted x and y arrays

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
    xs_arr = xs_arr[sort_idx]
    ys_arr = ys_arr[sort_idx]

    # If these are energies (rank 1 array), normalize to zero at far field
    if normalize_energy and ys_arr.ndim == 1:
        # shift to zero at largest separation (last after ascending sort)
        ys_arr = ys_arr - ys_arr[-1]

    return xs_arr, ys_arr


def calc_curve_diff_auc(
    seps_ref: ArrayLike,
    e_ref: ArrayLike,
    seps_pred: ArrayLike,
    e_pred: ArrayLike,
    *,
    seps_range: tuple[float | None, float | None] = (None, None),
    normalize: bool = True,
    interpolate: bool | int = False,
) -> float:
    """Calculate the absolute area under the curve of the difference between two curves.
    Handles different x-samplings by interpolating to a common grid.

    Args:
        seps_ref (ArrayLike[float]): Reference interatomic distances (Å)
        e_ref (ArrayLike[float]): Reference potential energies (eV)
        seps_pred (ArrayLike[float]): Predicted interatomic distances (Å)
        e_pred (ArrayLike[float]): Predicted potential energies (eV)
        seps_range (tuple[float | None, float | None] | None): Optional range of
            interatomic distances to consider. Can be None to auto-set based on data
            range. If tuple is None, uses intersection of both curves' x-ranges.
        normalize (bool): Whether to normalize by reference curve's bounding box area.
        interpolate (bool | int): If False (default), uses the provided points directly.
            If True, uses 100 points for interpolation.
            If an integer, uses that many points for interpolation.

    Returns:
        float: Absolute area under the curve of the difference between the curves.
            If normalize=True, returns unitless value, otherwise in eV·Å.
    """
    # Validate and sort both curves
    seps_ref, e_ref = _validate_diatomic_curve(seps_ref, e_ref, normalize_energy=False)
    seps_pred, e_pred = _validate_diatomic_curve(
        seps_pred, e_pred, normalize_energy=False
    )

    # Check if interpolation is needed
    if not interpolate and not np.array_equal(seps_ref, seps_pred):
        raise ValueError(
            f"Reference and predicted distances must be same when {interpolate=}\n"
            f"{seps_ref=}, {seps_pred=}"
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

    if interpolate:
        # Create grid for interpolation
        n_points = 100 if interpolate is True else interpolate
        seps_interp = np.linspace(seps_min, seps_max, n_points)

        # Interpolate both curves to the common grid
        e_ref_interp = np.interp(seps_interp, seps_ref, e_ref)
        e_pred_interp = np.interp(seps_interp, seps_pred, e_pred)

        # Calculate absolute difference and integrate
        diff = np.abs(e_ref_interp - e_pred_interp)
        auc = np.trapezoid(diff, seps_interp)
    else:
        # If no interpolation, restrict to the requested range
        mask = (seps_ref >= seps_min) & (seps_ref <= seps_max)
        if not np.any(mask):
            raise ValueError(f"No points within range {seps_min=}..{seps_max=}")
        seps_ref = seps_ref[mask]
        e_ref = e_ref[mask]
        e_pred = e_pred[mask]
        diff = np.abs(e_ref - e_pred)
        auc = np.trapezoid(diff, seps_ref)

    if normalize:
        # Normalize by bounding box of reference curve on the (possibly masked) domain.
        # When interpolate=True, uses full ref range; when False, uses masked subset.
        seps_span, e_span = np.ptp(seps_ref), np.ptp(e_ref)
        box_area = seps_span * e_span
        if box_area > 0:
            auc = auc / box_area

    # Ensure AUC is always positive
    return float(np.abs(auc))


def calc_energy_mae(
    seps_ref: ArrayLike,
    e_ref: ArrayLike,
    seps_pred: ArrayLike,
    e_pred: ArrayLike,
    *,
    interpolate: bool | int = False,
) -> float:
    """Calculate mean absolute error between two energy curves.
    Handles different x-samplings by interpolating to a common grid.

    Args:
        seps_ref (ArrayLike[float]): Reference interatomic distances (Å)
        e_ref (ArrayLike[float]): Reference potential energies (eV)
        seps_pred (ArrayLike[float]): Predicted interatomic distances (Å)
        e_pred (ArrayLike[float]): Predicted potential energies (eV)
        interpolate (bool | int): If False (default), uses the provided points directly.
            If True, uses 100 points for interpolation.
            If an integer, uses that many points for interpolation.

    Returns:
        float: Mean absolute error between the curves (eV).
    """
    # Validate and sort both curves
    seps_ref, e_ref = _validate_diatomic_curve(seps_ref, e_ref, normalize_energy=False)
    seps_pred, e_pred = _validate_diatomic_curve(
        seps_pred, e_pred, normalize_energy=False
    )

    # Check if interpolation is needed
    if not interpolate and not np.array_equal(seps_ref, seps_pred):
        raise ValueError(
            f"Reference and predicted distances must be same when {interpolate=}\n"
            f"{seps_ref=}, {seps_pred=}"
        )

    # Get data range bounds
    data_min = max(seps_ref.min(), seps_pred.min())
    data_max = min(seps_ref.max(), seps_pred.max())

    if interpolate:
        # Create grid for interpolation
        n_points = 100 if interpolate is True else interpolate
        seps_interp = np.linspace(data_min, data_max, n_points)

        # Interpolate both curves to the common grid
        e_ref_interp = np.interp(seps_interp, seps_ref, e_ref)
        e_pred_interp = np.interp(seps_interp, seps_pred, e_pred)

        # Calculate MAE on interpolated data
        return float(np.mean(np.abs(e_ref_interp - e_pred_interp)))
    # If no interpolation, calculate MAE directly
    return float(np.mean(np.abs(e_ref - e_pred)))


def calc_second_deriv_smoothness(seps: ArrayLike, energies: ArrayLike) -> float:
    """Calculate smoothness using RMS of second derivative (lower is smoother)."""
    seps_arr, energies_arr = _validate_diatomic_curve(
        seps, energies, normalize_energy=False
    )
    d2y = np.gradient(np.gradient(energies_arr, seps_arr), seps_arr)  # ty: ignore[no-matching-overload]
    return float(np.sqrt(np.mean(d2y**2)))


def calc_total_variation_smoothness(seps: ArrayLike, energies: ArrayLike) -> float:
    """Calculate smoothness using mean absolute gradient (lower is smoother)."""
    seps, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)
    dy = np.gradient(energies, seps)
    return float(np.log10(np.mean(np.abs(dy))))


def calc_curvature_smoothness(seps: ArrayLike, energies: ArrayLike) -> float:
    """Calculate smoothness using mean absolute curvature (lower is smoother)."""
    seps, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)
    dy = np.gradient(energies, seps)
    d2y = np.gradient(dy, seps)
    curvature = np.abs(d2y) / (1 + dy**2) ** 1.5
    return float(np.log10(np.mean(curvature)))


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
    _, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)

    tv_energy = np.sum(np.abs(np.diff(energies)))
    e_min = np.min(energies)
    direct_energy_diff = abs(energies[0] - e_min) + abs(energies[-1] - e_min)

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
    _, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)
    _, _, flips = _threshold_diff_signs(energies)
    return float(np.sum(flips))


def calc_energy_grad_norm_max(seps: ArrayLike, energies: ArrayLike) -> float:
    """Calculate maximum absolute value of energy gradient.

    Args:
        seps (ArrayLike[float]): Interatomic distances in Å.
        energies (ArrayLike[float]): Energies in eV.

    Returns:
        float: Maximum absolute value of energy gradient.
    """
    seps, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)
    return float(np.max(np.abs(np.gradient(energies, seps))))


def calc_energy_jump(seps: ArrayLike, energies: ArrayLike) -> float:
    """Calculate energy jump metric as sum of absolute energy differences at flip
    points.

    Args:
        seps (ArrayLike[float]): Interatomic distances in Å.
        energies (ArrayLike[float]): Energies in eV.

    Returns:
        float: Sum of absolute energy differences at flip points.
    """
    _, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)
    diffs, _, flips = _threshold_diff_signs(energies)
    return float(np.abs(diffs[:-1][flips]).sum() + np.abs(diffs[1:][flips]).sum())
