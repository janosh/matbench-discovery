"""Test force-based metrics for diatomic curves."""

from collections.abc import Callable

import numpy as np
import pytest

from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves


def x_forces(first_atom_x: list[float] | np.ndarray) -> np.ndarray:
    """Forces (n_distances, 2, 3) with the given x-component on the first atom and
    the equal and opposite force on the second.
    """
    forces = np.zeros((len(first_atom_x), 2, 3))
    forces[:, 0, 0] = first_atom_x
    forces[:, 1, 0] = -forces[:, 0, 0]
    return forces


OSCILLATING = [1.0, 2.0, -1.0, 3.0, -2.0]  # diffs [1, -3, 4, -5]
STEP = [1.0] * 20 + [-1.0] * 20


@pytest.mark.parametrize(
    ("metric_fn", "first_atom_x", "kwargs", "expected"),
    [
        (diatomics.calc_force_flips, OSCILLATING, {}, 3),
        (diatomics.calc_force_flips, STEP, {}, 1),
        # sub-threshold forces count as zero, so only a stricter threshold sees them
        (diatomics.calc_force_flips, [0.005, -0.005, 1, -1], {}, 1),
        (diatomics.calc_force_flips, [0.005, -0.005, 1, -1], {"threshold": 1e-3}, 3),
        # |1| + |-3| + |4| + |-5|
        (diatomics.calc_force_total_variation, OSCILLATING, {}, 13),
        (diatomics.calc_force_total_variation, np.linspace(-1, 1, 40), {}, 2),
        # all 3 diffs flip: (|1|+|3|) + (|3|+|4|) + (|4|+|5|)
        (diatomics.calc_force_jump, OSCILLATING, {}, 20),
        (diatomics.calc_force_jump, [1.0, 2.0, 3.0, 4.0, 5.0], {}, 0),
    ],
)
def test_force_smoothness_metrics(
    metric_fn: Callable[..., float],
    first_atom_x: list[float] | np.ndarray,
    kwargs: dict[str, float],
    expected: float,
) -> None:
    """Force flips, total variation and jump on hand-computed curves."""
    seps = np.linspace(5, 1, len(first_atom_x))
    assert metric_fn(seps, x_forces(first_atom_x), **kwargs) == pytest.approx(
        expected, abs=1e-12
    )


def test_force_mae(
    pred_ref_diatomic_curves: tuple[DiatomicCurves, DiatomicCurves],
) -> None:
    """Force MAE on a shared grid and interpolation across mismatched grids."""
    ref_curves, pred_curves = pred_ref_diatomic_curves
    ref_curve, pred_curve = ref_curves.homo_nuclear["H"], pred_curves.homo_nuclear["H"]
    x_ref, f_ref = ref_curve.distances, ref_curve.forces
    x_pred, f_pred = pred_curve.distances, pred_curve.forces

    mae = diatomics.calc_force_mae(x_ref, f_ref, x_pred, f_pred)
    # same 240 |diffs| summed in reverse (sorted) order: measured rel error 1.9e-16,
    # pairwise-summation bound ~log2(240) ~ 8 f64 eps
    assert mae == pytest.approx(np.mean(np.abs(f_ref - f_pred)), rel=1.8e-15, abs=0)

    x_pred_modified = x_pred * 1.05  # 5% different spacing
    with pytest.raises(
        ValueError, match="Reference and predicted distances must be same"
    ):
        diatomics.calc_force_mae(
            x_ref, f_ref, x_pred_modified, f_pred, interpolate=False
        )

    # default and custom interpolation grids sample the same overlap consistently
    mae_interp = diatomics.calc_force_mae(
        x_ref, f_ref, x_pred_modified, f_pred, interpolate=True
    )
    mae_custom_interp = diatomics.calc_force_mae(
        x_ref, f_ref, x_pred_modified, f_pred, interpolate=50
    )
    assert abs(mae_interp - mae_custom_interp) < 0.5


def test_force_mae_interpolation_uses_shared_distance_range() -> None:
    """Interpolation should compare force curves on their overlapping input domain."""
    x_ref = np.array([2.0, 3.0, 4.0])
    x_pred = np.array([2.1, 3.1, 4.1])

    mae = diatomics.calc_force_mae(
        x_ref, x_forces(x_ref), x_pred, x_forces(x_pred), interpolate=20
    )
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
