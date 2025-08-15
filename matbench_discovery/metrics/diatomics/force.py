"""Force-based metrics for diatomic curves."""

import numpy as np
from numpy.typing import ArrayLike

from matbench_discovery.metrics.diatomics.energy import _validate_diatomic_curve


def calc_force_mae(
    seps_ref: ArrayLike,
    f_ref: np.ndarray,
    seps_pred: ArrayLike,
    f_pred: np.ndarray,
    *,
    interpolate: bool | int = False,
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
        interpolate (bool | int): If False (default), uses the provided points directly.
            If True, uses 100 points for interpolation.
            If an integer, uses that many points for interpolation.

    Returns:
        float: Mean absolute error between the curves (eV/Å).
    """
    # Validate and sort both curves
    seps_ref, f_ref = _validate_diatomic_curve(seps_ref, f_ref, normalize_energy=False)
    seps_pred, f_pred = _validate_diatomic_curve(
        seps_pred, f_pred, normalize_energy=False
    )

    # Check if interpolation is needed
    if not interpolate:
        # If no interpolation is needed and distances match, calculate MAE directly
        if np.array_equal(seps_ref, seps_pred):
            return float(np.mean(np.abs(f_ref - f_pred)))
        raise ValueError(
            f"Reference and predicted distances must be same when {interpolate=}\n"
            f"{seps_ref=}, {seps_pred=}"
        )

    # Create grid for interpolation
    n_points = 100 if interpolate is True else interpolate
    seps_interp = np.logspace(1, -1, n_points)

    # Interpolate each component separately
    f_ref_interp = np.zeros((len(seps_interp), *f_ref.shape[1:]))
    f_pred_interp = np.zeros((len(seps_interp), *f_pred.shape[1:]))
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
    seps: ArrayLike,
    forces: np.ndarray,
    threshold: float = 1e-2,  # 10meV/A threshold as in reference code
) -> float:
    """Calculate number of (unphysical) force direction changes.

    Args:
        seps (Sequence[float]): Interatomic distances in Å.
        forces (np.ndarray): Forces of shape (n_distances, n_atoms, 3).
        threshold (float, optional): Forces below this threshold (in eV/Å) are
            considered zero. Defaults to 1e-2 (10 meV/Å).

    Returns:
        float: Number of force direction changes.
    """
    # Sort by separations in descending order
    _, forces = _validate_diatomic_curve(seps, forces, normalize_energy=False)

    fs = forces[:, 0, 0]

    rounded_fs = np.copy(fs)
    rounded_fs[np.abs(rounded_fs) < threshold] = 0
    fs_sign = np.sign(rounded_fs)
    mask = fs_sign != 0
    rounded_fs = rounded_fs[mask]
    fs_sign = fs_sign[mask]
    f_flip = np.diff(fs_sign) != 0
    return float(np.sum(f_flip))


def calc_force_total_variation(seps: ArrayLike, forces: np.ndarray) -> float:
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


def calc_force_jump(seps: ArrayLike, forces: np.ndarray) -> float:
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
    seps: ArrayLike,
    energies: ArrayLike,
    forces: np.ndarray,  # shape (n_distances, n_atoms, 3)
    *,
    interpolate: bool | int = False,
) -> float:
    """Calculate mean absolute deviation between forces and -dE/dr.

    Args:
        seps (ArrayLike): Interatomic distances in Å.
        energies (ArrayLike): Energies in eV.
        forces (np.ndarray): Forces acting on atoms at each separation of shape
            (n_distances, n_atoms, 3).
        interpolate (bool | int): If False (default), uses the provided points directly.
            If True, uses 100 points for interpolation.
            If an integer, uses that many points for interpolation.

    Returns:
        float: Mean absolute deviation between forces and -dE/dr.
    """
    seps, energies = _validate_diatomic_curve(seps, energies, normalize_energy=False)
    seps, forces = _validate_diatomic_curve(seps, forces, normalize_energy=False)

    if interpolate:
        # Create grid for interpolation
        n_points = 100 if interpolate is True else int(interpolate)
        seps_interp = np.linspace(seps.min(), seps.max(), n_points)

        # Interpolate energies
        energies_interp = np.interp(seps_interp, seps, energies)

        # Interpolate forces (only x-component for simplicity)
        forces_interp = np.zeros((len(seps_interp), forces.shape[1], forces.shape[2]))
        for atom_idx in range(forces.shape[1]):
            for dim in range(forces.shape[2]):
                forces_interp[:, atom_idx, dim] = np.interp(
                    seps_interp, seps, forces[:, atom_idx, dim]
                )

        # Calculate energy gradient using central differences on interpolated data
        energy_grad = np.gradient(energies_interp, seps_interp)
        # Compare forces with energy gradient
        return float(np.mean(np.abs(forces_interp + energy_grad.reshape(-1, 1, 1))))
    # Calculate energy gradient using central differences
    energy_grad = np.gradient(energies, seps)

    # Compare only x-component of forces with energy gradient
    # For diatomic molecules, forces should be equal and opposite
    # on the two atoms along the x-axis
    return float(np.mean(np.abs(forces + energy_grad.reshape(-1, 1, 1))))
