"""Force-based metrics for diatomic curves."""

from collections.abc import Sequence

import numpy as np

from matbench_discovery.metrics.diatomics.energy import _validate_diatomic_curve


def calc_force_mae(
    seps_ref: Sequence[float],
    f_ref: np.ndarray,
    seps_pred: Sequence[float],
    f_pred: np.ndarray,
) -> float:
    """Calculate mean absolute error between two force curves.
    Handles different x-samplings by interpolating to a common grid.

    Args:
        seps_ref (Sequence[float]): Reference interatomic distances (Å)
        f_ref (np.ndarray): Reference forces of shape
            (n_distances, n_atoms, 3)
        seps_pred (Sequence[float]): Predicted interatomic distances (Å)
        f_pred (np.ndarray): Predicted forces of shape
            (n_distances, n_atoms, 3)

    Returns:
        float: Mean absolute error between the curves (eV/Å).
    """
    # Validate and sort both curves
    seps_ref, f_ref = _validate_diatomic_curve(seps_ref, f_ref)
    seps_pred, f_pred = _validate_diatomic_curve(seps_pred, f_pred)

    # Get data range bounds
    data_min = max(seps_ref.min(), seps_pred.min())
    data_max = min(seps_ref.max(), seps_pred.max())

    # Create a fine grid for interpolation
    seps_interp = np.linspace(data_min, data_max, 1000)

    # Initialize interpolated arrays
    f_ref_interp = np.zeros((len(seps_interp), *f_ref.shape[1:]))
    f_pred_interp = np.zeros((len(seps_interp), *f_pred.shape[1:]))

    # Interpolate each component separately
    for atom_idx in range(f_ref.shape[1]):
        for dim in range(3):
            f_ref_interp[:, atom_idx, dim] = np.interp(
                seps_interp, seps_ref, f_ref[:, atom_idx, dim]
            )
            f_pred_interp[:, atom_idx, dim] = np.interp(
                seps_interp, seps_pred, f_pred[:, atom_idx, dim]
            )

    # Calculate MAE
    return float(np.mean(np.abs(f_ref_interp - f_pred_interp)))


def calc_force_flips(
    seps: Sequence[float],
    forces: np.ndarray,
    threshold: float = 1e-2,  # 10meV/A threshold as in reference code
) -> int:
    """Calculate number of (unphysical) force direction changes.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        forces (np.ndarray): Forces of shape (n_distances, n_atoms, 3).
        threshold (float, optional): Forces below this threshold (in eV/Å) are
            considered zero. Defaults to 1e-2 (10 meV/Å).

    Returns:
        int: Number of force direction changes.
    """
    seps = np.asarray(seps)
    # Use x-component of force on first atom
    forces_x = forces[:, 0, 0]

    # Round forces near zero (avoid numerical sensitivity)
    rounded_fs = np.copy(forces_x)
    rounded_fs[np.abs(rounded_fs) < threshold] = 0
    fs_sign = np.sign(rounded_fs)

    # Mask out zero values
    mask = fs_sign != 0
    fs_sign = fs_sign[mask]

    # Count sign changes
    return int(np.sum(np.diff(fs_sign) != 0))


def calc_force_total_variation(seps: Sequence[float], forces: np.ndarray) -> float:
    """Calculate total variation in forces.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        forces (np.ndarray): Forces of shape (n_distances, n_atoms, 3).

    Returns:
        float: Sum of absolute differences between consecutive force values.
    """
    sort_idx = np.argsort(seps)[::-1]  # sort in descending order
    forces_x = forces[sort_idx, 0, 0]  # x-component of force on first atom
    return float(np.sum(np.abs(np.diff(forces_x))))


def calc_force_jump(seps: Sequence[float], forces: np.ndarray) -> float:
    """Calculate force jump metric as sum of absolute force differences at flip points.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        forces (np.ndarray): Forces of shape (n_distances, n_atoms, 3).

    Returns:
        float: Sum of absolute force differences at flip points.
    """
    sort_idx = np.argsort(seps)[::-1]  # sort in descending order
    forces_x = forces[sort_idx, 0, 0]  # x-component of force on first atom

    f_diff = np.diff(forces_x)
    f_diff_sign = np.sign(f_diff)
    mask = f_diff_sign != 0
    f_diff = f_diff[mask]
    f_diff_sign = f_diff_sign[mask]
    f_diff_flip = np.diff(f_diff_sign) != 0

    force_jumps = np.abs(f_diff[:-1][f_diff_flip]).sum() + np.abs(
        f_diff[1:][f_diff_flip]
    )
    return float(force_jumps.sum())


def calc_conservation_deviation(
    seps: Sequence[float],
    energies: Sequence[float],
    forces: np.ndarray,  # shape (n_distances, n_atoms, 3)
) -> float:
    """Calculate mean absolute deviation between forces and -dE/dr.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        energies (Sequence[float]): Energies in eV.
        forces (np.ndarray): Forces acting on atoms at each separation of shape
            (n_distances, n_atoms, 3).

    Returns:
        float: Mean absolute deviation between forces and -dE/dr.
    """
    _sorted_seps, energies = _validate_diatomic_curve(seps, energies)
    seps, forces = _validate_diatomic_curve(seps, forces)

    # Calculate energy gradient using central differences
    energy_grad = np.gradient(energies, seps)

    # Compare only x-component of forces with energy gradient
    # For diatomic molecules, forces should be equal and opposite
    # on the two atoms along the x-axis
    return float(np.mean(np.abs(forces + energy_grad.reshape(-1, 1, 1))))
