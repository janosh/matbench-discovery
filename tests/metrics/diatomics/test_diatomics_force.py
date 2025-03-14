"""Test force-based metrics for diatomic curves."""

import re

import numpy as np
import pytest

from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves

PredRefForces = tuple[
    dict[str, tuple[np.ndarray, np.ndarray]], dict[str, tuple[np.ndarray, np.ndarray]]
]


@pytest.fixture
def pred_ref_forces(
    pred_ref_diatomic_curves: tuple[DiatomicCurves, DiatomicCurves],
) -> PredRefForces:
    """Create test data for diatomic curves.

    Returns:
        PredRefForces: Reference and predicted curves for each element.
    """
    ref_curves, pred_curves = pred_ref_diatomic_curves
    ref_homo_nuc, pred_homo_nuc = ref_curves.homo_nuclear, pred_curves.homo_nuclear
    f_h_ref, f_h_pred = ref_homo_nuc["H"].forces, pred_homo_nuc["H"].forces
    f_he_ref, f_he_pred = ref_homo_nuc["He"].forces, pred_homo_nuc["He"].forces
    ref_dists, pred_dists = ref_curves.distances, pred_curves.distances

    ref_forces = {"H": (ref_dists, f_h_ref), "He": (ref_dists, f_he_ref)}
    pred_forces = {"H": (pred_dists, f_h_pred), "He": (pred_dists, f_he_pred)}
    return ref_forces, pred_forces


def test_force_flips(pred_ref_forces: PredRefForces) -> None:
    """Test force flips calculation."""
    _, pred_curves = pred_ref_forces
    dists, forces = pred_curves["H"]

    # Test with default parameters
    flips = diatomics.calc_force_flips(dists, forces)
    assert flips == 20  # Updated for logspace distances

    # Test with different force threshold
    flips_strict = diatomics.calc_force_flips(dists, forces, threshold=1e-3)
    assert flips_strict >= flips  # Stricter threshold should find more flips

    # Test with simple cases
    xs = np.logspace(1, -1, 40)
    # Create forces array with shape (n_distances, n_atoms, 3)
    forces = np.zeros((len(xs), 2, 3))
    # Force on first atom changes sign once
    forces[:20, 0, 0] = 1
    forces[20:, 0, 0] = -1
    # Force on second atom is equal and opposite
    forces[:, 1, 0] = -forces[:, 0, 0]

    flips = diatomics.calc_force_flips(xs, forces)
    assert flips == 1  # One sign change


def test_force_total_variation(pred_ref_forces: PredRefForces) -> None:
    """Test force total variation calculation."""
    _, pred_curves = pred_ref_forces
    dists, forces = pred_curves["H"]

    # Test with default parameters
    total_var = diatomics.calc_force_total_variation(dists, forces)
    assert isinstance(total_var, float)
    assert total_var >= 0  # Total variation should be non-negative

    # Test with simple cases
    dists = np.logspace(1, -1, 40)
    # Create forces array with shape (n_distances, n_atoms, 3)
    forces = np.zeros((len(dists), 2, 3))
    # Linear force on first atom
    forces[:, 0, 0] = np.linspace(-1, 1, len(dists))
    # Force on second atom is equal and opposite
    forces[:, 1, 0] = -forces[:, 0, 0]

    total_var = diatomics.calc_force_total_variation(dists, forces)
    # Total variation of linear function from -1 to 1
    assert total_var == pytest.approx(2.0, abs=0.1)

    # For logspace distances, the gradient of a linear energy function
    # will not be constant, so we can't expect total variation to be 0
    dists = np.logspace(1, -1, 40)
    e_linear = np.linspace(-1, 1, len(dists))
    f_linear = np.gradient(e_linear, dists).reshape(-1, 1, 1)
    # Just check that it's a reasonable value
    assert diatomics.calc_force_total_variation(dists, f_linear) == pytest.approx(
        4.0, abs=0.5
    )


def test_force_jump(pred_ref_forces: PredRefForces) -> None:
    """Test force jump calculation."""
    _, pred_curves = pred_ref_forces
    dists, forces = pred_curves["H"]

    # Test with default parameters
    jump = diatomics.calc_force_jump(dists, forces)
    assert isinstance(jump, float)
    assert jump >= 0  # Force jump should be non-negative


def test_force_mae(pred_ref_forces: PredRefForces) -> None:
    """Test force MAE calculation."""
    ref_forces, pred_forces = pred_ref_forces
    x_ref, f_ref = ref_forces["H"]
    x_pred, f_pred = pred_forces["H"]

    # Test with default parameters
    mae = diatomics.calc_force_mae(x_ref, f_ref, x_pred, f_pred)
    assert isinstance(mae, float)
    assert mae >= 0  # MAE should be non-negative


def test_force_mae_interpolation(pred_ref_forces: PredRefForces) -> None:
    """Test interpolation behavior in force MAE calculation."""
    ref_forces, pred_forces = pred_ref_forces
    x_ref, f_ref = ref_forces["H"]
    x_pred, f_pred = pred_forces["H"]

    # Create modified x_pred with different spacing
    x_pred_modified = x_pred * 1.05  # 5% difference

    # Test with interpolation=False (should raise error when x values don't match)
    with pytest.raises(
        ValueError, match="Reference and predicted distances must be same"
    ):
        diatomics.calc_force_mae(
            x_ref, f_ref, x_pred_modified, f_pred, interpolate=False
        )

    # Test with interpolation=True
    mae_interp = diatomics.calc_force_mae(
        x_ref, f_ref, x_pred_modified, f_pred, interpolate=True
    )
    assert isinstance(mae_interp, float)
    assert mae_interp >= 0  # MAE should be non-negative

    # Test with custom number of interpolation points
    mae_custom_interp = diatomics.calc_force_mae(
        x_ref, f_ref, x_pred_modified, f_pred, interpolate=50
    )
    assert isinstance(mae_custom_interp, float)
    assert mae_custom_interp >= 0  # MAE should be non-negative

    # results should be similar but not identical due to different interpolation grids
    assert mae_interp != mae_custom_interp
    assert abs(mae_interp - mae_custom_interp) < 0.5  # Should be reasonably close


def test_invalid_inputs() -> None:
    """Test that invalid inputs raise appropriate errors."""
    seps = np.array([1, 2])
    energies = np.array([0, 1])
    # Forces shape: (n_seps, n_atoms=2, xyz=3)
    forces = np.zeros((2, 2, 3))

    # Test with NaN values
    with pytest.raises(ValueError, match="Input contains NaN values"):
        diatomics.calc_conservation_deviation(np.array([np.nan, 2]), energies, forces)

    # Test with infinite values
    with pytest.raises(ValueError, match="Input contains infinite values"):
        diatomics.calc_conservation_deviation(np.array([np.inf, 2]), energies, forces)

    # Test with duplicate x values
    with pytest.raises(ValueError, match="xs contains 1 duplicates"):
        diatomics.calc_conservation_deviation(np.array([1, 1]), energies, forces)

    # Test with mismatched array sizes
    with pytest.raises(ValueError, match=re.escape("len(xs)=2 != len(ys)=1")):
        diatomics.calc_conservation_deviation(seps, np.array([0]), forces)

    # Test with single point
    with pytest.raises(ValueError, match="Input must have at least 2 points"):
        diatomics.calc_conservation_deviation(
            np.array([1]), np.array([0]), np.zeros((1, 2, 3))
        )


def test_conservation_deviation(
    pred_ref_diatomic_curves: tuple[DiatomicCurves, DiatomicCurves],
) -> None:
    """Test conservation deviation calculation."""
    # Simple parabolic potential
    seps = np.logspace(1, -1, 40)
    energies = 0.5 * (seps - 2) ** 2
    # Forces should be F = -dE/dr = -(r-2)
    # Shape: (n_seps, n_atoms=2, xyz=3)
    forces = np.zeros((len(seps), 2, 3))
    # Set x-component of forces for both atoms (equal and opposite)
    forces[:, 0, 0] = -(seps - 2)  # force on first atom
    forces[:, 1, 0] = seps - 2  # force on second atom

    deviation = diatomics.calc_conservation_deviation(seps, energies, forces)
    assert deviation == pytest.approx(2.0, abs=0.5)

    _, pred_curves = pred_ref_diatomic_curves
    h_data = pred_curves.homo_nuclear["H"]
    dists, e_pred, _f_pred = h_data.distances, h_data.energies, h_data.forces

    # Test with default parameters
    # Create forces array with correct shape
    forces_pred = np.zeros((len(dists), 2, 3))
    grad = np.gradient(e_pred, dists)
    forces_pred[:, 0, 0] = -grad  # force on first atom
    forces_pred[:, 1, 0] = grad  # force on second atom

    dev = diatomics.calc_conservation_deviation(dists, e_pred, forces_pred)
    assert isinstance(dev, float)
    assert dev >= 0  # Deviation should be non-negative

    # Test with sorted vs unsorted input (should give same result)
    np_rng = np.random.default_rng(seed=0)
    perm = np_rng.permutation(len(dists))
    dists_unsorted = dists[perm]
    e_unsorted = e_pred[perm]
    forces_unsorted = forces_pred[perm]
    dev_unsorted = diatomics.calc_conservation_deviation(
        dists_unsorted, e_unsorted, forces_unsorted
    )
    assert dev == pytest.approx(dev_unsorted)

    dists = np.arange(1, 6)
    e_linear = np.arange(1, 6)
    # Constant curve
    e_const = np.ones(len(dists))

    # Create forces arrays with correct shape
    forces_linear = np.zeros((len(dists), 2, 3))
    forces_linear[:, 0, 0] = -1  # constant force on first atom
    forces_linear[:, 1, 0] = 1  # constant force on second atom

    forces_const = np.zeros((len(dists), 2, 3))  # zero forces everywhere

    linear_cons = diatomics.calc_conservation_deviation(dists, e_linear, forces_linear)
    const_cons = diatomics.calc_conservation_deviation(dists, e_const, forces_const)
    # Linear curve should have constant force equal to slope=1
    assert linear_cons == pytest.approx(1, abs=1e-3)
    # Constant curve should have zero force
    assert const_cons == pytest.approx(0, abs=1e-3)

    seps = np.logspace(1, -1, 40)
    energies = np.zeros_like(seps)  # flat potential
    forces = np.zeros((len(seps), 2, 3))  # zero forces

    # Test conservation deviation
    assert diatomics.calc_conservation_deviation(
        seps, energies, forces
    ) == pytest.approx(0, abs=1e-3)


def test_conservation_deviation_interpolation(pred_ref_forces: PredRefForces) -> None:
    """Test interpolation behavior in conservation deviation calculation."""
    _, pred_curves = pred_ref_forces
    dists, forces = pred_curves["H"]

    # Create energies from forces by integration
    energies = np.zeros_like(dists)
    for i in range(1, len(dists)):
        # Simple trapezoidal integration
        dx = dists[i - 1] - dists[i]
        f_avg = (forces[i - 1, 0, 0] + forces[i, 0, 0]) / 2
        energies[i] = energies[i - 1] - f_avg * dx

    # Test with interpolation=False (default)
    deviation = diatomics.calc_conservation_deviation(
        dists, energies, forces, interpolate=False
    )
    assert isinstance(deviation, float)
    assert deviation >= 0  # Deviation should be non-negative

    # Test with interpolation=True
    deviation_interp = diatomics.calc_conservation_deviation(
        dists, energies, forces, interpolate=True
    )
    assert isinstance(deviation_interp, float)
    assert deviation_interp >= 0  # Deviation should be non-negative

    # Test with custom number of interpolation points
    deviation_custom_interp = diatomics.calc_conservation_deviation(
        dists, energies, forces, interpolate=50
    )
    assert isinstance(deviation_custom_interp, float)
    assert deviation_custom_interp >= 0  # Deviation should be non-negative

    # results should be similar but not identical due to different interpolation grids
    assert deviation_interp != deviation_custom_interp
    assert (
        abs(deviation_interp - deviation_custom_interp) < 0.5
    )  # Should be reasonably close
