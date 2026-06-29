"""Test diatomic curve metrics calculation functions."""

import re
from collections.abc import Callable

import numpy as np
import pytest

from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves
from matbench_discovery.metrics.diatomics.energy import (
    calc_energy_diff_flips,
    calc_energy_grad_norm_max,
    calc_energy_jump,
    calc_energy_mae,
)

EnergyCurve = tuple[np.ndarray, np.ndarray]
PredRefEnergies = tuple[EnergyCurve, EnergyCurve]


@pytest.fixture
def pred_ref_e_curves(
    pred_ref_diatomic_curves: tuple[DiatomicCurves, DiatomicCurves],
) -> PredRefEnergies:
    """Return reference and predicted H energy curves from the shared fixture."""
    ref_curves, pred_curves = pred_ref_diatomic_curves
    ref_curve = ref_curves.distances, ref_curves.homo_nuclear["H"].energies
    pred_curve = pred_curves.distances, pred_curves.homo_nuclear["H"].energies
    return ref_curve, pred_curve


@pytest.mark.parametrize(
    "metric_func, kwargs, max_result",
    [
        (diatomics.calc_curve_diff_auc, {}, 10.0),
        (diatomics.calc_curve_diff_auc, {"seps_range": (1.0, 3.0)}, 10.0),
        (calc_energy_mae, {}, None),
    ],
    ids=["auc-default", "auc-range", "energy-mae"],
)
def test_curve_pair_metrics_on_fixture(
    metric_func: Callable[..., float],
    kwargs: dict[str, object],
    max_result: float | None,
    pred_ref_e_curves: PredRefEnergies,
) -> None:
    """Curve-pair energy metrics return finite non-negative floats on fixture data."""
    (seps_ref, energies_ref), (seps_pred, energies_pred) = pred_ref_e_curves

    result = metric_func(seps_ref, energies_ref, seps_pred, energies_pred, **kwargs)
    assert isinstance(result, float)
    assert result >= 0
    if max_result is not None:
        assert result < max_result


@pytest.mark.parametrize(
    "metric_func, custom_interpolate, max_delta, max_result",
    [
        (diatomics.calc_curve_diff_auc, 500, 0.5, 10.0),
        (calc_energy_mae, 110, 1.0, None),
    ],
)
def test_curve_pair_metric_interpolation(
    metric_func: Callable[..., float],
    custom_interpolate: int,
    max_delta: float,
    max_result: float | None,
    pred_ref_e_curves: PredRefEnergies,
) -> None:
    """Curve-pair metrics handle interpolation consistently."""
    (seps_ref, energies_ref), (seps_pred, energies_pred) = pred_ref_e_curves

    # Create modified x_pred with different spacing
    seps_pred_modified = seps_pred * 1.05  # 5% difference

    # Test with interpolation=False (should raise error when x values don't match)
    with pytest.raises(
        ValueError,
        match="Reference and predicted distances must be same",
    ):
        metric_func(
            seps_ref,
            energies_ref,
            seps_pred_modified,
            energies_pred,
            interpolate=False,
        )

    # Test with interpolation=True (default)
    metric_interp = metric_func(
        seps_ref, energies_ref, seps_pred_modified, energies_pred, interpolate=True
    )
    assert isinstance(metric_interp, float)
    assert metric_interp >= 0
    if max_result is not None:
        assert metric_interp < max_result

    # Test with custom number of interpolation points
    metric_custom_interp = metric_func(
        seps_ref,
        energies_ref,
        seps_pred_modified,
        energies_pred,
        interpolate=custom_interpolate,
    )
    assert isinstance(metric_custom_interp, float)
    assert metric_custom_interp >= 0
    if max_result is not None:
        assert metric_custom_interp < max_result

    # results should be similar but not identical due to different interpolation grids
    assert metric_interp != metric_custom_interp
    assert abs(metric_interp - metric_custom_interp) < max_delta


@pytest.mark.parametrize(
    "metric_func",
    [
        diatomics.calc_tortuosity,
        diatomics.calc_energy_jump,
        diatomics.calc_energy_diff_flips,
        diatomics.calc_energy_grad_norm_max,
    ],
)
def test_energy_metrics_on_fixture(
    metric_func: Callable[..., float],
    pred_ref_e_curves: PredRefEnergies,
) -> None:
    """All single-curve energy metrics return non-negative floats on fixture data."""
    _, pred_curve = pred_ref_e_curves
    seps_pred, energies_pred = pred_curve
    result = metric_func(seps_pred, energies_pred)
    assert isinstance(result, float)
    assert result >= 0


@pytest.mark.parametrize(
    "vals, threshold, expected_n_flips, expected_diffs_len",
    [
        # Monotone increasing: no flips
        (np.array([1.0, 2.0, 3.0, 4.0]), 1e-3, 0, 3),
        # Single oscillation: 1 flip (up then down)
        (np.array([1.0, 3.0, 2.0]), 1e-3, 1, 2),
        # Multiple oscillations: 3 flips
        (np.array([1.0, 3.0, 2.0, 4.0, 1.0]), 1e-3, 3, 4),
        # Diffs below threshold get zeroed → no flips
        (np.array([1.0, 1.0004, 1.0002, 1.0006]), 1e-3, 0, 0),
        # Two points: 1 diff, 0 flips
        (np.array([1.0, 2.0]), 1e-3, 0, 1),
    ],
)
def test_threshold_diff_signs(
    vals: np.ndarray, threshold: float, expected_n_flips: int, expected_diffs_len: int
) -> None:
    """Test _threshold_diff_signs with known inputs."""
    from matbench_discovery.metrics.diatomics.energy import _threshold_diff_signs

    diffs, signs, flips = _threshold_diff_signs(vals, threshold)
    assert len(diffs) == expected_diffs_len
    assert len(signs) == expected_diffs_len
    assert int(np.sum(flips)) == expected_n_flips
    # signs should only be -1 or +1 (zeros are filtered)
    if len(signs) > 0:
        assert set(signs.tolist()) <= {-1.0, 1.0}


@pytest.mark.parametrize(
    "energies, expected_flips, expected_jump",
    [
        # Monotone: 0 flips, 0 jump
        (np.array([1.0, 2.0, 3.0, 4.0, 5.0]), 0, 0.0),
        # diffs=[2,-1,2,1], signs=[+,-,+,+], flips at positions 0,1 → 2 flips
        # jump = |diffs[0]|+|diffs[1]| + |diffs[1]|+|diffs[2]| = 2+1+1+2 = 6
        (np.array([1.0, 3.0, 2.0, 4.0, 5.0]), 2, 6.0),
        # diffs=[2,-1,2,-2.5], signs=[+,-,+,-], 3 flips
        # jump = |2|+|1| + |1|+|2| + |2|+|2.5| = 3+3+4.5 = 10.5
        (np.array([0.0, 2.0, 1.0, 3.0, 0.5]), 3, 10.5),
    ],
)
def test_energy_flips_and_jumps_concrete(
    energies: np.ndarray, expected_flips: int, expected_jump: float
) -> None:
    """Test calc_energy_diff_flips and calc_energy_jump with concrete values."""
    seps = np.arange(1, len(energies) + 1, dtype=float)
    assert calc_energy_diff_flips(seps, energies) == expected_flips
    assert calc_energy_jump(seps, energies) == pytest.approx(expected_jump)


@pytest.mark.parametrize(
    "energies, expected_grad_max",
    [
        (np.array([5.0, 5.0, 5.0, 5.0]), 0.0),  # constant: gradient = 0
        # Linear slope=2: gradient = 2 everywhere
        (np.array([2.0, 4.0, 6.0, 8.0]), 2.0),
        # Quadratic x^2 at x=1..4: max gradient at endpoint = 7
        (np.array([1.0, 4.0, 9.0, 16.0]), 7.0),
    ],
)
def test_energy_grad_norm_max_concrete(
    energies: np.ndarray, expected_grad_max: float
) -> None:
    """Test calc_energy_grad_norm_max with analytically known gradients."""
    seps = np.arange(1, len(energies) + 1, dtype=float)
    assert calc_energy_grad_norm_max(seps, energies) == pytest.approx(expected_grad_max)


@pytest.mark.parametrize(
    "seps, energies_ref, energies_pred, normalize, expected",
    [
        (
            np.array([1.0, 2.0, 3.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            True,
            0.0,
        ),
        (
            np.array([1.0, 2.0, 3.0, 4.0]),
            np.zeros(4),
            np.ones(4),
            False,
            3.0,
        ),
        (
            np.array([1.0, 2.0, 3.0, 4.0]),
            np.array([0.0, 2.0, 0.0, 2.0]),  # ptp=2, so box_area = 3 * 2 = 6
            np.zeros(4),
            True,
            0.5,
        ),
    ],
    ids=["identical", "raw", "normalized"],
)
def test_curve_diff_auc_concrete_cases(
    seps: np.ndarray,
    energies_ref: np.ndarray,
    energies_pred: np.ndarray,
    normalize: bool,
    expected: float,
) -> None:
    """Test AUC exact values for identical, raw, and normalized cases."""
    auc = diatomics.calc_curve_diff_auc(
        seps, energies_ref, seps, energies_pred, normalize=normalize
    )
    assert auc == pytest.approx(expected)


def test_validate_normalize_energy() -> None:
    """Test that normalize_energy shifts energies to zero at far field."""
    from matbench_discovery.metrics.diatomics.energy import _validate_diatomic_curve

    seps = np.array([1.0, 2.0, 3.0, 4.0])
    energies = np.array([10.0, 5.0, 2.0, 1.0])
    _, normed = _validate_diatomic_curve(seps, energies, normalize_energy=True)
    # After sort ascending, last value (at sep=4) should be 0
    assert normed[-1] == 0.0
    assert normed[0] == pytest.approx(9.0)  # 10 - 1

    # Forces (ndim > 1) should NOT be normalized even with flag set
    forces = np.ones((4, 2, 3))
    _, forces_out = _validate_diatomic_curve(seps, forces, normalize_energy=True)
    np.testing.assert_array_equal(forces_out, forces)


@pytest.mark.parametrize(
    "energies, expected",
    [
        (np.arange(1, 6), 1.0),
        (np.arange(1, 6) ** 2, 1.0),
        (np.ones(5), np.nan),
        (np.array([1, 1, 1.5, 2, 2]), 1.0),
    ],
    ids=["linear", "quadratic", "constant", "step"],
)
def test_tortuosity_simple_cases(energies: np.ndarray, expected: float) -> None:
    """Test tortuosity with simple, easy-to-reason-about cases."""
    dists = np.arange(1, 6)
    result = diatomics.calc_tortuosity(dists, energies)
    if np.isnan(expected):
        assert np.isnan(result)
    else:
        assert result == pytest.approx(expected)


@pytest.mark.parametrize(
    "xs, ys, match",
    [
        (np.array([1, np.nan, 3, 4, 5]), np.arange(1, 6), "Input contains NaN"),
        (np.arange(1, 6), np.array([1, np.inf, 3, 4, 5]), "Input contains infinite"),
        (np.array([1]), np.array([1]), "Input must have at least 2 points"),
        (np.array([1, 1, 2, 3, 4]), np.arange(1, 6), "xs contains 1 duplicates"),
        # non-adjacent duplicates (regression: np.diff used to miss these)
        (np.array([1, 3, 2, 1, 4]), np.arange(1, 6), "xs contains 1 duplicates"),
        (np.array([1, 2, 1, 3, 2]), np.arange(1, 6), "xs contains 2 duplicates"),
        (
            np.array([1, 2, 3]),
            np.arange(1, 6),
            re.escape("len(xs_arr)=3 != len(ys_arr)=5"),
        ),
    ],
)
def test_validation_errors(xs: np.ndarray, ys: np.ndarray, match: str) -> None:
    """Test that _validate_diatomic_curve rejects invalid inputs."""
    with pytest.raises(ValueError, match=match):
        calc_energy_grad_norm_max(xs, ys)


def test_unsorted_input_accepted() -> None:
    """Unsorted x values should work since sorting is handled internally."""
    xs = np.array([1, 3, 2, 5, 4])
    ys = np.array([1, 3, 2, 5, 4])
    assert np.isfinite(
        diatomics.calc_curve_diff_auc(np.arange(1, 6), np.arange(1, 6), xs, ys)
    )
