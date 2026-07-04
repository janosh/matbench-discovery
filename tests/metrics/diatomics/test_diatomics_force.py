"""Test force-based metrics for diatomic curves."""

import numpy as np
import pytest

from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves

PredRefForces = tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]


@pytest.fixture
def pred_ref_forces(
    pred_ref_diatomic_curves: tuple[DiatomicCurves, DiatomicCurves],
) -> PredRefForces:
    """Create test force curves for H-H."""
    ref_curves, pred_curves = pred_ref_diatomic_curves
    ref_curve = ref_curves.homo_nuclear["H"]
    pred_curve = pred_curves.homo_nuclear["H"]
    return (
        (ref_curve.distances, ref_curve.forces),
        (pred_curve.distances, pred_curve.forces),
    )


def test_force_flips(pred_ref_forces: PredRefForces) -> None:
    """Test force flips calculation."""
    _, (dists, forces) = pred_ref_forces

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
    _, (dists, forces) = pred_ref_forces

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
    _, (dists, forces) = pred_ref_forces

    jump = diatomics.calc_force_jump(dists, forces)
    assert isinstance(jump, float)
    assert jump >= 0

    # Concrete case: forces with one sign change → measurable jump
    seps = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    forces_osc = np.zeros((5, 2, 3))
    forces_osc[:, 0, 0] = [1.0, 2.0, -1.0, 3.0, -2.0]  # two flips
    jump_concrete = diatomics.calc_force_jump(seps, forces_osc)
    assert jump_concrete > 0
    # diffs=[1,-3,4,-5], signs=[+,-,+,-], all 3 flips
    # jump = |1|+|3| + |3|+|4| + |4|+|5| = 4+7+9 = 20
    assert jump_concrete == pytest.approx(20.0)

    # Monotone forces: zero jump
    forces_mono = np.zeros((5, 2, 3))
    forces_mono[:, 0, 0] = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert diatomics.calc_force_jump(seps, forces_mono) == pytest.approx(0.0)


def test_force_mae(pred_ref_forces: PredRefForces) -> None:
    """Force MAE with default parameters and interpolation across mismatched grids."""
    (x_ref, f_ref), (x_pred, f_pred) = pred_ref_forces

    # Test with default parameters
    mae = diatomics.calc_force_mae(x_ref, f_ref, x_pred, f_pred)
    assert isinstance(mae, float)
    assert mae >= 0  # MAE should be non-negative

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


def test_force_mae_interpolation_uses_shared_distance_range() -> None:
    """Interpolation should compare force curves on their overlapping input domain."""
    x_ref = np.array([2.0, 3.0, 4.0])
    x_pred = np.array([2.1, 3.1, 4.1])
    f_ref = np.zeros((3, 2, 3))
    f_pred = np.zeros((3, 2, 3))
    f_ref[:, 0, 0] = x_ref
    f_pred[:, 0, 0] = x_pred

    mae = diatomics.calc_force_mae(x_ref, f_ref, x_pred, f_pred, interpolate=20)
    assert mae == pytest.approx(0)


@pytest.mark.parametrize(
    "x_pred",
    [np.array([4.0, 5.0, 6.0]), np.array([3.0, 4.0, 5.0])],
    ids=["disjoint", "single_shared_point"],
)
def test_force_mae_interpolation_rejects_unusable_overlap(x_pred: np.ndarray) -> None:
    """Force MAE interpolation rejects disjoint or single-point-overlap curves."""
    x_ref = np.array([1.0, 2.0, 3.0])
    forces = np.zeros((3, 2, 3))

    with pytest.raises(ValueError, match="no overlap"):
        diatomics.calc_force_mae(x_ref, forces, x_pred, forces, interpolate=True)
