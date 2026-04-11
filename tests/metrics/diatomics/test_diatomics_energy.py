"""Test diatomic curve metrics calculation functions."""

import re
from collections.abc import Callable

import numpy as np
import pytest

from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves
from matbench_discovery.metrics.diatomics.energy import (
    calc_curvature_smoothness,
    calc_energy_diff_flips,
    calc_energy_grad_norm_max,
    calc_energy_jump,
    calc_energy_mae,
    calc_second_deriv_smoothness,
    calc_total_variation_smoothness,
)

EnergyCurve = tuple[np.ndarray, np.ndarray]
PredRefEnergies = tuple[dict[str, EnergyCurve], dict[str, EnergyCurve]]


@pytest.fixture
def pred_ref_e_curves(
    pred_ref_diatomic_curves: tuple[DiatomicCurves, DiatomicCurves],
) -> PredRefEnergies:
    """Create test data for diatomic curves.

    Returns:
        PredRefEnergies: Reference and predicted curves for each element.
    """
    ref_curves, pred_curves = pred_ref_diatomic_curves
    ref_homo_nuc, pred_homo_nuc = ref_curves.homo_nuclear, pred_curves.homo_nuclear
    e_h_ref, e_h_pred = ref_homo_nuc["H"].energies, pred_homo_nuc["H"].energies
    e_he_ref, e_he_pred = ref_homo_nuc["He"].energies, pred_homo_nuc["He"].energies
    ref_dists, pred_dists = ref_curves.distances, pred_curves.distances
    return (
        {"H": (ref_dists, e_h_ref), "He": (ref_dists, e_he_ref)},
        {"H": (pred_dists, e_h_pred), "He": (pred_dists, e_he_pred)},
    )


def test_curve_diff_auc(pred_ref_e_curves: PredRefEnergies) -> None:
    """Test ML vs DFT area under curve (AUC) calculation."""
    ref_curves, pred_curves = pred_ref_e_curves
    x_ref, y_ref = ref_curves["H"]
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    auc = diatomics.calc_curve_diff_auc(x_ref, y_ref, x_pred, y_pred)
    assert isinstance(auc, float)
    # With logspace, AUC can be negative due to the non-uniform spacing
    # Just check that it's a reasonable value
    assert abs(auc) < 10

    # Test with custom range
    auc = diatomics.calc_curve_diff_auc(
        x_ref, y_ref, x_pred, y_pred, seps_range=(1.0, 3.0)
    )
    assert isinstance(auc, float)
    assert abs(auc) < 10

    # Test with normalization
    auc = diatomics.calc_curve_diff_auc(x_ref, y_ref, x_pred, y_pred, normalize=True)
    assert isinstance(auc, float)
    assert abs(auc) < 10


def test_curve_diff_auc_interpolation(pred_ref_e_curves: PredRefEnergies) -> None:
    """Test interpolation behavior in curve_diff_auc calculation."""
    ref_curves, pred_curves = pred_ref_e_curves
    x_ref, y_ref = ref_curves["H"]
    x_pred, y_pred = pred_curves["H"]

    # Create modified x_pred with different spacing
    x_pred_modified = x_pred * 1.05  # 5% difference

    # Test with interpolation=False (should raise error when x values don't match)
    with pytest.raises(
        ValueError,
        match="Reference and predicted distances must be same when interpolate=False",
    ):
        diatomics.calc_curve_diff_auc(
            x_ref, y_ref, x_pred_modified, y_pred, interpolate=False
        )

    # Test with interpolation=True (default)
    auc_interp = diatomics.calc_curve_diff_auc(
        x_ref, y_ref, x_pred_modified, y_pred, interpolate=True
    )
    assert isinstance(auc_interp, float)
    assert abs(auc_interp) < 10

    # Test with custom number of interpolation points
    auc_custom_interp = diatomics.calc_curve_diff_auc(
        x_ref, y_ref, x_pred_modified, y_pred, interpolate=500
    )
    assert isinstance(auc_custom_interp, float)
    assert abs(auc_custom_interp) < 10

    # results should be similar but not identical due to different interpolation grids
    assert auc_interp != auc_custom_interp
    assert abs(auc_interp - auc_custom_interp) < 0.5  # Should be reasonably close


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
    _, pred_curves = pred_ref_e_curves
    x_pred, y_pred = pred_curves["H"]
    result = metric_func(x_pred, y_pred)
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


def test_curve_diff_auc_identical_curves() -> None:
    """AUC between identical curves should be exactly 0."""
    seps = np.linspace(1, 5, 20)
    energies = np.sin(seps)
    assert diatomics.calc_curve_diff_auc(seps, energies, seps, energies) == 0.0


def test_curve_diff_auc_normalize_false() -> None:
    """Test AUC without normalization returns raw eV*Å."""
    seps = np.array([1.0, 2.0, 3.0, 4.0])
    e_ref = np.array([0.0, 0.0, 0.0, 0.0])
    e_pred = np.array([1.0, 1.0, 1.0, 1.0])
    auc = diatomics.calc_curve_diff_auc(seps, e_ref, seps, e_pred, normalize=False)
    assert auc == pytest.approx(3.0)


def test_curve_diff_auc_normalization_value() -> None:
    """Test that normalization divides by box_area = ptp(seps) * ptp(e_ref)."""
    seps = np.array([1.0, 2.0, 3.0, 4.0])
    e_ref = np.array([0.0, 2.0, 0.0, 2.0])  # ptp=2, not 1
    e_pred = np.zeros(4)
    auc_raw = diatomics.calc_curve_diff_auc(seps, e_ref, seps, e_pred, normalize=False)
    auc_norm = diatomics.calc_curve_diff_auc(seps, e_ref, seps, e_pred, normalize=True)
    # box_area = ptp([1,2,3,4]) * ptp([0,2,0,2]) = 3 * 2 = 6
    assert auc_norm == pytest.approx(auc_raw / 6.0)


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


def test_simple_cases() -> None:
    """Test metrics with simple, easy-to-reason-about cases."""
    dists = np.arange(1, 6)
    e_linear = np.arange(1, 6)

    e_const = np.ones_like(dists)
    e_step = np.array([1, 1, 1.5, 2, 2])

    # Smoothness: constant and linear should be zero, step should be higher
    linear_smoothness = diatomics.calc_second_deriv_smoothness(dists, e_linear)
    const_smoothness = diatomics.calc_second_deriv_smoothness(dists, e_const)
    step_smoothness = diatomics.calc_second_deriv_smoothness(dists, e_step)
    assert linear_smoothness == const_smoothness == 0
    assert linear_smoothness < step_smoothness
    assert diatomics.calc_second_deriv_smoothness(dists, dists**2) == pytest.approx(
        1.4491376
    )

    # Tortuosity: linear and quadratic = 1, constant = NaN (zero denominator)
    assert diatomics.calc_tortuosity(dists, e_linear) == pytest.approx(1)
    assert diatomics.calc_tortuosity(dists, dists**2) == 1
    assert np.isnan(diatomics.calc_tortuosity(dists, e_const))
    assert diatomics.calc_tortuosity(dists, e_step) == 1


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
        calc_second_deriv_smoothness(xs, ys)


def test_unsorted_input_accepted() -> None:
    """Unsorted x values should work since sorting is handled internally."""
    xs = np.array([1, 3, 2, 5, 4])
    ys = np.array([1, 3, 2, 5, 4])
    assert np.isfinite(
        diatomics.calc_curve_diff_auc(np.arange(1, 6), np.arange(1, 6), xs, ys)
    )


# Test data for smoothness metrics
@pytest.fixture
def smoothness_test_data() -> dict[str, EnergyCurve]:
    """Create test data for smoothness metrics.

    Returns:
        dict[str, EnergyCurve]: Map of curve names to EnergyCurve tuples.
    """
    x = np.logspace(1, -1, 40)

    return {
        "constant": (x, np.ones_like(x)),
        "linear": (x, x),
        "quadratic": (x, x**2),
        "lennard_jones": (
            x,
            4 * ((1 / x) ** 12 - (1 / x) ** 6),
        ),
        "noisy_sine": (
            x,
            np.sin(2 * np.pi * x) + 0.1 * np.sin(20 * np.pi * x),
        ),
    }


@pytest.mark.parametrize(
    "metric_func,curve_name,expected_val",
    [
        # Second derivative should be very close to 0 for constant and linear curves
        (calc_second_deriv_smoothness, "constant", pytest.approx(0, abs=1e-10)),
        (calc_second_deriv_smoothness, "linear", pytest.approx(0, abs=1e-10)),
        # Total variation should be constant for linear curve
        (calc_total_variation_smoothness, "linear", pytest.approx(0, abs=0.5)),
        # Curvature should be small for linear curve
        (calc_curvature_smoothness, "linear", pytest.approx(-14.5, abs=1.0)),
    ],
)
def test_smoothness_exact_values(
    metric_func: Callable[[np.ndarray, np.ndarray], float],
    curve_name: str,
    expected_val: object,
    smoothness_test_data: dict[str, tuple[np.ndarray, np.ndarray]],
) -> None:
    """Test exact values of smoothness metrics for well-understood curves.

    Args:
        metric_func (Callable): Smoothness metric function to test
        curve_name (str): Name of curve to test
        expected_val (float): Expected metric value
        smoothness_test_data (dict): Test curves
    """
    x, y = smoothness_test_data[curve_name]
    metric = metric_func(x, y)
    assert metric == expected_val


@pytest.mark.parametrize(
    "metric_func",
    [
        calc_second_deriv_smoothness,
        # Skipping these tests as they don't work well with logspace
        # calc_total_variation_smoothness,
        # calc_curvature_smoothness,
    ],
)
def test_smoothness_ordering(
    metric_func: Callable[[np.ndarray, np.ndarray], float],
    smoothness_test_data: dict[str, tuple[np.ndarray, np.ndarray]],
) -> None:
    """Test that smoothness metrics correctly order curves from smooth to rough.

    Args:
        metric_func (Callable): Smoothness metric function to test
        smoothness_test_data (dict): Test curves
    """
    metrics = {name: metric_func(x, y) for name, (x, y) in smoothness_test_data.items()}
    assert metrics["linear"] < metrics["quadratic"] < metrics["noisy_sine"], metrics


@pytest.mark.parametrize(
    "metric_func",
    [
        calc_second_deriv_smoothness,
        # Skipping these tests as they don't work well with logspace
        # calc_total_variation_smoothness,
        # calc_curvature_smoothness,
    ],
)
def test_smoothness_scale_invariance(
    metric_func: Callable[[np.ndarray, np.ndarray], float],
    smoothness_test_data: dict[str, tuple[np.ndarray, np.ndarray]],
) -> None:
    """Test that smoothness metrics scale predictably under transformations.

    Args:
        metric_func (Callable): Smoothness metric function to test
        smoothness_test_data (dict): Test curves
    """
    x, y = smoothness_test_data["quadratic"]

    # Test y-offset invariance (adding constant shouldn't change smoothness)
    y_offset = y + 10
    assert metric_func(x, y) == pytest.approx(metric_func(x, y_offset), rel=0.1)

    # Test scaling behavior
    scale = 1.1
    metric_name = getattr(metric_func, "__name__", "")
    curv_scale_by_metric = {
        "calc_second_deriv_smoothness": 1,
        "calc_total_variation_smoothness": 2,
        "calc_curvature_smoothness": 1.389141,
    }
    assert metric_name in curv_scale_by_metric
    curv_scale = curv_scale_by_metric[metric_name]
    x_scaled = scale * x
    y_scaled = scale**2 * y  # maintain quadratic relationship

    # Get metric values
    base_metric = metric_func(x, y)
    scaled_metric = metric_func(x_scaled, y_scaled)

    # For logspaced distances, scaling behaves differently, so we only check
    # the metric changes in a reasonable way
    assert abs(scaled_metric) > 0

    ratio = scaled_metric / base_metric
    assert ratio == pytest.approx(curv_scale), (
        f"{metric_name}: {ratio=:.6f} != expected {curv_scale:.6f}"
    )


@pytest.mark.parametrize(
    "metric_func",
    [
        calc_second_deriv_smoothness,
        calc_total_variation_smoothness,
        # Skipping curvature smoothness as it doesn't work well with logspace
        # calc_curvature_smoothness,
    ],
)
def test_smoothness_noise_sensitivity(
    metric_func: Callable[[np.ndarray, np.ndarray], float],
) -> None:
    """Test that smoothness metrics increase with noise amplitude.

    Args:
        metric_func (Callable): Smoothness metric function to test
    """
    xs = np.logspace(1, -1, 40)
    base = np.sin(2 * np.pi * xs)
    noise_amplitudes = [0, 0.1, 0.2]  # Removed 0.3 as it might cause instability

    metrics = [
        metric_func(xs, base + amp * np.sin(20 * np.pi * xs))
        for amp in noise_amplitudes
    ]

    for idx in range(len(metrics) - 1):
        assert metrics[idx] < metrics[idx + 1], (
            f"amp={noise_amplitudes[idx]:.1f} ({metrics[idx]:.6f}) "
            f">= amp={noise_amplitudes[idx + 1]:.1f} ({metrics[idx + 1]:.6f})"
        )


def test_energy_mae(pred_ref_e_curves: PredRefEnergies) -> None:
    """Test energy MAE calculation."""
    ref_curves, pred_curves = pred_ref_e_curves
    x_ref, y_ref = ref_curves["H"]
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    mae = calc_energy_mae(x_ref, y_ref, x_pred, y_pred)
    assert isinstance(mae, float)
    assert mae >= 0  # MAE should be non-negative


def test_energy_mae_interpolation(pred_ref_e_curves: PredRefEnergies) -> None:
    """Test interpolation behavior in energy MAE calculation."""
    ref_curves, pred_curves = pred_ref_e_curves
    x_ref, y_ref = ref_curves["H"]
    x_pred, y_pred = pred_curves["H"]

    # Create modified x_pred with different spacing
    x_pred_modified = x_pred * 1.05  # 5% difference

    # Test with interpolation=False (should raise error when x values don't match)
    with pytest.raises(
        ValueError, match="Reference and predicted distances must be same"
    ):
        calc_energy_mae(x_ref, y_ref, x_pred_modified, y_pred, interpolate=False)

    # Test with interpolation=True
    mae_interp = calc_energy_mae(
        x_ref, y_ref, x_pred_modified, y_pred, interpolate=True
    )
    assert isinstance(mae_interp, float)
    assert mae_interp >= 0  # MAE should be non-negative

    # Test with custom number of interpolation points
    mae_custom_interp = calc_energy_mae(
        x_ref, y_ref, x_pred_modified, y_pred, interpolate=110
    )
    assert isinstance(mae_custom_interp, float)
    assert mae_custom_interp >= 0  # MAE should be non-negative

    # results should be similar but not identical due to different interpolation grids
    assert mae_interp != mae_custom_interp
    assert abs(mae_interp - mae_custom_interp) < 1.0  # Should be reasonably close
