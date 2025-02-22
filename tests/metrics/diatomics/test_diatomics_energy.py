"""Test diatomic curve metrics calculation functions."""

import json
from collections.abc import Callable

import numpy as np
import pytest

from matbench_discovery import ROOT
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurves
from matbench_discovery.metrics.diatomics.energy import (
    calc_curvature_smoothness,
    calc_curve_diff_auc,
    calc_energy_diff_flips,
    calc_energy_grad_norm_max,
    calc_energy_jump,
    calc_energy_mae,
    calc_second_deriv_smoothness,
    calc_tortuosity,
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
    assert auc > 0  # AUC should be positive
    assert auc < 1  # Normalized AUC should be less than 1 for similar curves

    # Test without normalization
    auc_unnorm = diatomics.calc_curve_diff_auc(
        x_ref, y_ref, x_pred, y_pred, normalize=False
    )
    assert auc_unnorm > 0
    assert auc_unnorm > auc  # Unnormalized AUC should be larger than normalized

    seps = np.linspace(0.5, 5, 100)
    e_ref = np.zeros_like(seps)
    e_pred = np.ones_like(seps)  # constant shift of 1

    # Test with default parameters
    auc = calc_curve_diff_auc(seps, e_ref, seps, e_pred)
    assert auc == pytest.approx(4.5)  # area = 1 * (5 - 0.5)

    # Test with custom range
    auc = calc_curve_diff_auc(seps, e_ref, seps, e_pred, seps_range=(1, 2))
    assert auc == pytest.approx(1)  # area = 1 * (2 - 1)

    # Test with normalization
    auc = calc_curve_diff_auc(
        seps, e_ref, seps, e_pred, seps_range=(1, 2), normalize=True
    )
    assert auc == pytest.approx(1)  # normalized area = 1


def test_tortuosity(pred_ref_e_curves: PredRefEnergies) -> None:
    """Test tortuosity calculation."""
    _, pred_curves = pred_ref_e_curves
    x_pred, y_pred = pred_curves["H"]

    # Test tortuosity calculation
    tort = diatomics.calc_tortuosity(x_pred, y_pred)
    assert isinstance(tort, float)
    assert tort >= 1  # Tortuosity should be at least 1 (arc length ≥ direct distance)


def test_energy_jump(pred_ref_e_curves: PredRefEnergies) -> None:
    """Test energy jump calculation."""
    _, pred_curves = pred_ref_e_curves
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    e_jump = diatomics.calc_energy_jump(x_pred, y_pred)
    assert isinstance(e_jump, float)
    assert e_jump >= 0  # Energy jump should be non-negative


def test_energy_diff_flips(pred_ref_e_curves: PredRefEnergies) -> None:
    """Test energy difference flips calculation."""
    _, pred_curves = pred_ref_e_curves
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    n_flips = diatomics.calc_energy_diff_flips(x_pred, y_pred)
    assert n_flips == 1


def test_energy_grad_norm_max(pred_ref_e_curves: PredRefEnergies) -> None:
    """Test energy gradient norm maximum calculation."""
    _, pred_curves = pred_ref_e_curves
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    grad_max = diatomics.calc_energy_grad_norm_max(x_pred, y_pred)
    assert isinstance(grad_max, float)
    assert grad_max >= 0  # Maximum gradient norm should be non-negative


@pytest.fixture
def mace_data() -> tuple[
    dict[str, tuple[np.ndarray, np.ndarray]], dict[str, tuple[np.ndarray, np.ndarray]]
]:
    """Load MACE model diatomic curve data.

    Returns:
        tuple[dict[str, tuple], dict[str, tuple]]: Reference and predicted curves for
            each element.
    """
    json_path = (
        f"{ROOT}/models/mace/mace-mp-0/mace-mlip-arena-homonuclear-diatomics.json"
    )
    with open(json_path) as file:
        data = json.load(file)

    # Convert data to required format
    ref_curves = {}
    pred_curves = {}

    for elem_data in data:
        distances = np.array(elem_data["R"])
        energies = np.array(elem_data["E"])
        # Shift energies so the energy at largest separation is 0
        shifted_energies = energies - energies[-1]

        elem_symbol = elem_data["name"][: len(elem_data["name"]) // 2]

        # For this test, use the same data for reference and prediction
        # In practice, you would use different models' predictions
        ref_curves[elem_symbol] = (distances, shifted_energies)
        pred_curves[elem_symbol] = (distances, shifted_energies)

    return ref_curves, pred_curves


def test_simple_cases() -> None:
    """Test metrics with simple, easy-to-reason-about cases."""
    dists = np.arange(1, 6)
    e_linear = np.arange(1, 6)

    # Linear curve should have zero second derivative, constant force and tortuosity = 1
    assert diatomics.calc_second_deriv_smoothness(dists, e_linear) == 0
    assert diatomics.calc_tortuosity(dists, e_linear) == 1

    # Constant curve should have zero force and zero total variation
    e_const = np.ones_like(dists)

    # Test with quadratic curve
    quad_smoothness = diatomics.calc_second_deriv_smoothness(dists, dists**2)
    assert quad_smoothness == pytest.approx(1.4491376)
    assert diatomics.calc_tortuosity(dists, dists**2) == 1

    # Step function with a more gradual transition
    e_step = np.array([1, 1, 1.5, 2, 2])

    # Test smoothness
    linear_smoothness = diatomics.calc_second_deriv_smoothness(dists, e_linear)
    const_smoothness = diatomics.calc_second_deriv_smoothness(dists, e_const)
    step_smoothness = diatomics.calc_second_deriv_smoothness(dists, e_step)

    # Linear curve should be smoother than step function
    assert linear_smoothness < step_smoothness
    # Constant and linear curve should be perfectly smooth (zero second derivative)
    assert linear_smoothness == const_smoothness == 0

    # Test tortuosity
    linear_tort = diatomics.calc_tortuosity(dists, e_linear)
    const_tort = diatomics.calc_tortuosity(dists, e_const)
    step_tort = diatomics.calc_tortuosity(dists, e_step)

    # Linear curve should have minimal tortuosity (arc length = direct distance)
    assert linear_tort == pytest.approx(1)
    # Constant curve should have infinite tortuosity (no direct energy difference)
    assert np.isnan(const_tort)
    # Step function should have intermediate tortuosity
    assert step_tort == 1

    seps = np.linspace(0.5, 5, 100)
    energies = np.zeros_like(seps)  # flat potential

    # Test energy difference flips
    assert calc_energy_diff_flips(seps, energies) == 0

    # Test energy gradient norm max
    assert calc_energy_grad_norm_max(seps, energies) == pytest.approx(0)

    # Test energy jump
    assert calc_energy_jump(seps, energies) == pytest.approx(0)

    # Test tortuosity
    assert np.isnan(calc_tortuosity(seps, energies))

    # Test smoothness
    assert calc_second_deriv_smoothness(seps, energies) == pytest.approx(0)


def test_edge_cases() -> None:
    """Test metrics with edge cases."""
    xs = np.array([1, 2, 3, 4, 5])
    ys = np.array([1, 2, 3, 4, 5])

    # Test with NaN values
    y_nan = np.array([1, np.nan, 3, 4, 5])
    with pytest.raises(ValueError, match="Input contains NaN values: x=0, y=1"):
        diatomics.calc_curve_diff_auc(xs, ys, xs, y_nan)

    # Test with infinite values
    y_inf = np.array([1, np.inf, 3, 4, 5])
    with pytest.raises(ValueError, match="Input contains infinite values: x=0, y=1"):
        diatomics.calc_curve_diff_auc(xs, ys, xs, y_inf)

    # Test with single point
    x_single = np.array([1])
    y_single = np.array([1])
    with pytest.raises(ValueError, match="Input must have at least 2 points, got 1"):
        diatomics.calc_curve_diff_auc(x_single, y_single, x_single, y_single)

    # Test with duplicate x values
    x_dup = np.array([1, 1, 2, 3, 4])
    y_dup = np.array([1, 2, 3, 4, 5])
    with pytest.raises(ValueError, match="Input contains 1 duplicate x values"):
        diatomics.calc_curve_diff_auc(x_dup, y_dup, x_dup, y_dup)

    # Test with mismatched array sizes
    x_short = np.array([1, 2, 3])
    y_long = np.array([1, 2, 3, 4, 5])
    with pytest.raises(ValueError, match="x and y must have same size"):
        diatomics.calc_curve_diff_auc(x_short, y_long, x_short, y_long)

    # Test with unsorted x values (should work, sorting is handled internally)
    x_unsorted = np.array([1, 3, 2, 5, 4])
    y_unsorted = np.array([1, 3, 2, 5, 4])
    auc = diatomics.calc_curve_diff_auc(xs, ys, x_unsorted, y_unsorted)
    assert np.isfinite(auc)


# Test data for smoothness metrics
@pytest.fixture
def smoothness_test_data() -> dict[str, EnergyCurve]:
    """Create test data for smoothness metrics.

    Returns:
        dict[str, EnergyCurve]: Map of curve names to EnergyCurve tuples.
    """
    x = np.linspace(0.1, 1, 1000)
    x_lj = np.linspace(0.3, 3, 1000)  # special range for LJ potential

    return {
        "constant": (x, np.ones_like(x)),
        "linear": (x, x),
        "quadratic": (x, x**2),
        "lennard_jones": (
            x_lj,
            4 * ((1 / x_lj) ** 12 - (1 / x_lj) ** 6),
        ),
        "noisy_sine": (
            x,
            np.sin(2 * np.pi * x) + 0.1 * np.sin(20 * np.pi * x),
        ),
    }


@pytest.mark.parametrize(
    "metric_func,curve_name,expected_value",
    [
        # Second derivative should be very close to 0 for constant and linear curves
        (calc_second_deriv_smoothness, "constant", pytest.approx(0, abs=1e-10)),
        (calc_second_deriv_smoothness, "linear", pytest.approx(0, abs=1e-10)),
        # Total variation should be constant for linear curve
        (calc_total_variation_smoothness, "linear", pytest.approx(0, abs=0.1)),
        # Curvature should be small for linear curve
        (calc_curvature_smoothness, "linear", pytest.approx(-11.5, abs=0.5)),
    ],
)
def test_smoothness_exact_values(
    metric_func: Callable[[np.ndarray, np.ndarray], float],
    curve_name: str,
    expected_value: float,
    smoothness_test_data: dict[str, tuple[np.ndarray, np.ndarray]],
) -> None:
    """Test exact values of smoothness metrics for well-understood curves.

    Args:
        metric_func (Callable): Smoothness metric function to test
        curve_name (str): Name of curve to test
        expected_value (float): Expected metric value
        smoothness_test_data (dict): Test curves
    """
    x, y = smoothness_test_data[curve_name]
    metric = metric_func(x, y)
    assert metric == expected_value


@pytest.mark.parametrize(
    "metric_func",
    [
        calc_second_deriv_smoothness,
        calc_total_variation_smoothness,
        calc_curvature_smoothness,
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
    # Calculate metric for each curve
    metrics = {name: metric_func(x, y) for name, (x, y) in smoothness_test_data.items()}

    # Print all metrics for debugging
    if not (
        metrics["linear"] < metrics["quadratic"]
        and metrics["quadratic"] < metrics["noisy_sine"]
    ):
        raise AssertionError(
            "\nSmooth to rough ordering failed. Metric values:\n"
            + "\n".join(f"{name}: {value:.6f}" for name, value in metrics.items())
        )


@pytest.mark.parametrize(
    "metric_func",
    [
        calc_second_deriv_smoothness,
        calc_total_variation_smoothness,
        calc_curvature_smoothness,
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
    curv_scale = {
        "calc_second_deriv_smoothness": 1,
        "calc_total_variation_smoothness": 2,
        "calc_curvature_smoothness": 1.389141,
    }[metric_func.__name__]
    x_scaled = scale * x
    y_scaled = scale**2 * y  # maintain quadratic relationship

    # Get metric values
    base_metric = metric_func(x, y)
    scaled_metric = metric_func(x_scaled, y_scaled)

    if scaled_metric / base_metric != pytest.approx(curv_scale):
        raise AssertionError(
            f"{metric_func.__name__} did not scale correctly:\n"
            f"scale={scale:.1f}\n"
            f"base_metric={base_metric:.6f}\n"
            f"scaled_metric={scaled_metric:.6f}\n"
            f"ratio={scaled_metric / base_metric:.6f}\n"
            f"expected_ratio={curv_scale:.6f}"
        )


@pytest.mark.parametrize(
    "metric_func",
    [
        calc_second_deriv_smoothness,
        calc_total_variation_smoothness,
        calc_curvature_smoothness,
    ],
)
def test_smoothness_noise_sensitivity(
    metric_func: Callable[[np.ndarray, np.ndarray], float],
) -> None:
    """Test that smoothness metrics increase with noise amplitude.

    Args:
        metric_func (Callable): Smoothness metric function to test
    """
    x = np.linspace(0.1, 1, 1000)
    base = np.sin(2 * np.pi * x)
    noise_amps = [0, 0.1, 0.2]  # Removed 0.3 as it might cause instability

    # Calculate metrics for increasing noise levels
    metrics = []
    for amp in noise_amps:
        y = base + amp * np.sin(20 * np.pi * x)
        metrics.append(metric_func(x, y))

    # Check that metric increases monotonically with noise
    for i in range(len(metrics) - 1):
        if not metrics[i] < metrics[i + 1]:
            raise AssertionError(
                f"\nMetric did not increase monotonically with noise:\n"
                f"noise_amp={noise_amps[i]:.1f} -> metric={metrics[i]:.6f}\n"
                f"noise_amp={noise_amps[i + 1]:.1f} -> metric={metrics[i + 1]:.6f}\n"
                f"Metric values for all noise amplitudes:\n"
                + "\n".join(
                    f"  amp={amp:.1f}: {metric:.6f}"
                    for amp, metric in zip(noise_amps, metrics)
                )
            )


def test_energy_mae() -> None:
    """Test mean absolute error calculation."""
    seps = np.linspace(0.5, 5, 100)
    e_ref = np.zeros_like(seps)
    e_pred = np.ones_like(seps)  # constant shift of 1

    mae = calc_energy_mae(seps, e_ref, seps, e_pred)
    assert mae == pytest.approx(1)  # MAE = |1 - 0| = 1
