"""Force-based metrics for diatomic curves."""

import numpy as np
from numpy.typing import ArrayLike

from matbench_discovery.metrics.diatomics.energy import (
    _common_grid_pair,
    _threshold_diff_signs,
    _validate_diatomic_curve,
)


def calc_force_mae(
    seps_ref: ArrayLike,
    f_ref: ArrayLike,
    seps_pred: ArrayLike,
    f_pred: ArrayLike,
    *,
    interpolate: bool | int = False,
) -> float:
    """Calculate mean absolute error between two force curves.
    Handles different x-samplings by interpolating to a common grid.

    Args:
        seps_ref (ArrayLike): Reference interatomic distances (Å)
        f_ref (ArrayLike): Reference forces of shape
            (n_distances, n_atoms, 3)
        seps_pred (ArrayLike): Predicted interatomic distances (Å)
        f_pred (ArrayLike): Predicted forces of shape
            (n_distances, n_atoms, 3)
        interpolate (bool | int): If False (default), uses the provided points directly.
            If True, uses 100 points for interpolation.
            If an integer, uses that many points for interpolation.

    Returns:
        float: Mean absolute error between the curves (eV/Å).
    """
    _, f_ref, f_pred = _common_grid_pair(
        seps_ref, f_ref, seps_pred, f_pred, interpolate=interpolate
    )
    return float(np.mean(np.abs(f_ref - f_pred)))


def calc_force_flips(
    seps: ArrayLike,
    forces: np.ndarray,
    threshold: float = 1e-2,  # 10meV/A threshold as in reference code
) -> float:
    """Calculate number of (unphysical) force direction changes.

    Args:
        seps (ArrayLike): Interatomic distances in Å.
        forces (np.ndarray): Forces of shape (n_distances, n_atoms, 3).
        threshold (float, optional): Forces below this threshold (in eV/Å) are
            considered zero. Defaults to 1e-2 (10 meV/Å).

    Returns:
        float: Number of force direction changes.
    """
    _, forces = _validate_diatomic_curve(seps, forces)

    fs = forces[:, 0, 0].copy()
    fs[np.abs(fs) < threshold] = 0
    fs_sign = np.sign(fs[fs != 0])
    return float(np.sum(np.diff(fs_sign) != 0))


def calc_force_total_variation(seps: ArrayLike, forces: np.ndarray) -> float:
    """Calculate total variation in forces.

    Args:
        seps (ArrayLike): Interatomic distances in Å.
        forces (np.ndarray): Forces of shape (n_distances, n_atoms, 3).

    Returns:
        float: Sum of absolute differences between consecutive force values.
    """
    _, forces = _validate_diatomic_curve(seps, forces)
    forces_x = forces[:, 0, 0]  # x-component of force on first atom
    return float(np.sum(np.abs(np.diff(forces_x))))


def calc_force_jump(seps: ArrayLike, forces: np.ndarray) -> float:
    """Calculate force jump metric as sum of absolute force differences at flip points.

    Args:
        seps (ArrayLike): Interatomic distances in Å.
        forces (np.ndarray): Forces of shape (n_distances, n_atoms, 3).

    Returns:
        float: Sum of absolute force differences at flip points.
    """
    _, forces = _validate_diatomic_curve(seps, forces)
    diffs, _, flips = _threshold_diff_signs(forces[:, 0, 0], threshold=0)
    return float(np.abs(diffs[:-1][flips]).sum() + np.abs(diffs[1:][flips]).sum())
