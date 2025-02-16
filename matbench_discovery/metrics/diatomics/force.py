"""Force-based metrics for diatomic curves."""

from collections.abc import Sequence

import numpy as np

from matbench_discovery.metrics.diatomics.energy import _validate_curve_input


def calc_force_mae_vs_ref(
    seps_ref: Sequence[float],
    f_ref: Sequence[float],
    seps_pred: Sequence[float],
    f_pred: Sequence[float],
) -> float:
    """Calculate mean absolute error between two force curves.
    Handles different x-samplings by interpolating to a common grid.

    Args:
        seps_ref (Sequence[float]): Reference interatomic distances (Å)
        f_ref (Sequence[float]): Reference forces (eV/Å)
        seps_pred (Sequence[float]): Predicted interatomic distances (Å)
        f_pred (Sequence[float]): Predicted forces (eV/Å)

    Returns:
        float: Mean absolute error between the curves (eV/Å).
    """
    # Validate and sort both curves
    seps_ref, f_ref = _validate_curve_input(seps_ref, f_ref)
    seps_pred, f_pred = _validate_curve_input(seps_pred, f_pred)

    # Get data range bounds
    data_min = max(seps_ref.min(), seps_pred.min())
    data_max = min(seps_ref.max(), seps_pred.max())

    # Create a fine grid for interpolation
    seps_interp = np.linspace(data_min, data_max, 1000)

    # Interpolate both curves to the common grid
    f_ref_interp = np.interp(seps_interp, seps_ref, f_ref)
    f_pred_interp = np.interp(seps_interp, seps_pred, f_pred)

    # Calculate MAE
    return float(np.mean(np.abs(f_ref_interp - f_pred_interp)))


def calc_force_flips(
    seps: Sequence[float],
    forces: Sequence[float],
    threshold: float = 1e-2,  # 10meV/A threshold as in reference code
) -> float:
    """Calculate number of unphysical force direction changes.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        forces (Sequence[float]): Forces in eV/Å.
        threshold (float, optional): Forces below this threshold (in eV/Å) are
            considered zero. Defaults to 1e-2 (10 meV/Å).

    Returns:
        float: Number of force direction changes.
    """
    seps, forces = map(np.asarray, (seps, forces))

    # Round forces near zero (avoid numerical sensitivity)
    rounded_fs = np.copy(forces)
    rounded_fs[np.abs(rounded_fs) < threshold] = 0
    fs_sign = np.sign(rounded_fs)

    # Mask out zero values
    mask = fs_sign != 0
    fs_sign = fs_sign[mask]

    # Count sign changes
    return float(np.sum(np.diff(fs_sign) != 0))


def calc_force_total_variation(seps: Sequence[float], forces: Sequence[float]) -> float:
    """Calculate total variation in forces.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        forces (Sequence[float]): Forces in eV/Å.

    Returns:
        float: Sum of absolute differences between consecutive force values.
    """
    seps, forces = map(np.asarray, (seps, forces))
    sort_idx = np.argsort(seps)[::-1]  # sort in descending order
    forces = forces[sort_idx]
    return float(np.sum(np.abs(np.diff(forces))))


def calc_force_jump(seps: Sequence[float], forces: Sequence[float]) -> float:
    """Calculate force jump metric as sum of absolute force differences at flip points.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        forces (Sequence[float]): Forces in eV/Å.

    Returns:
        float: Sum of absolute force differences at flip points.
    """
    seps, forces = map(np.asarray, (seps, forces))
    sort_idx = np.argsort(seps)[::-1]  # sort in descending order
    forces = forces[sort_idx]

    fdiff = np.diff(forces)
    fdiff_sign = np.sign(fdiff)
    mask = fdiff_sign != 0
    fdiff = fdiff[mask]
    fdiff_sign = fdiff_sign[mask]
    fdiff_flip = np.diff(fdiff_sign) != 0

    return float(
        np.abs(fdiff[:-1][fdiff_flip]).sum() + np.abs(fdiff[1:][fdiff_flip]).sum()
    )
