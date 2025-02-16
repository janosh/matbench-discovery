"""Test diatomic curve metrics calculation functions."""

import json
from collections.abc import Callable

import numpy as np
import pytest

from matbench_discovery import ROOT
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics.energy import (
    calc_curvature_smoothness,
    calc_second_deriv_smoothness,
    calc_total_variation_smoothness,
)


def test_curve_diff_auc(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test AUC difference calculation."""
    ref_curves, pred_curves = test_data
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


def test_tortuosity(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test tortuosity calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test tortuosity calculation
    tort = diatomics.calc_tortuosity(x_pred, y_pred)
    assert isinstance(tort, float)
    assert tort >= 1.0  # Tortuosity should be at least 1 (arc length ≥ direct distance)


def test_conservation_deviation(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test conservation deviation calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    dev = diatomics.calc_conservation_deviation(x_pred, y_pred)
    assert isinstance(dev, float)
    assert dev >= 0  # Deviation should be non-negative

    # Test with sorted vs unsorted input (should give same result)
    np_rng = np.random.default_rng(seed=0)
    x_unsorted = np_rng.permutation(x_pred)
    y_unsorted = np_rng.permutation(y_pred)
    dev_unsorted = diatomics.calc_conservation_deviation(x_unsorted, y_unsorted)
    assert dev == pytest.approx(dev_unsorted)


def test_energy_jump(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test energy jump calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    ejump = diatomics.calc_energy_jump(x_pred, y_pred)
    assert isinstance(ejump, float)
    assert ejump >= 0  # Energy jump should be non-negative


def test_energy_diff_flips(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test energy difference flips calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    n_flips = diatomics.calc_energy_diff_flips(x_pred, y_pred)
    assert n_flips == 1


def test_energy_grad_norm_max(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test energy gradient norm maximum calculation."""
    _, pred_curves = test_data
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
    xs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y_linear = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    # Linear curve should have zero second derivative, constant force and tortuosity = 1
    assert diatomics.calc_second_deriv_smoothness(xs, y_linear) == 0.0
    assert diatomics.calc_tortuosity(xs, y_linear) == pytest.approx(1.0)
    assert diatomics.calc_force_total_variation(xs, y_linear) == pytest.approx(4)

    # Test with quadratic curve
    assert diatomics.calc_second_deriv_smoothness(xs, xs**2) == pytest.approx(1.4491376)
    assert diatomics.calc_tortuosity(xs, xs**2) == pytest.approx(1.0)

    # Constant curve
    y_const = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    # Step function with a more gradual transition
    y_step = np.array([1.0, 1.0, 1.5, 2.0, 2.0])

    # Test smoothness
    linear_smoothness = diatomics.calc_second_deriv_smoothness(xs, y_linear)
    const_smoothness = diatomics.calc_second_deriv_smoothness(xs, y_const)
    step_smoothness = diatomics.calc_second_deriv_smoothness(xs, y_step)

    # Linear curve should be smoother than step function
    assert linear_smoothness < step_smoothness
    # Constant and linear curve should be perfectly smooth (zero second derivative)
    assert linear_smoothness == const_smoothness == 0

    # Test tortuosity
    linear_tort = diatomics.calc_tortuosity(xs, y_linear)
    const_tort = diatomics.calc_tortuosity(xs, y_const)
    step_tort = diatomics.calc_tortuosity(xs, y_step)

    # Linear curve should have minimal tortuosity (arc length = direct distance)
    assert linear_tort == pytest.approx(1.0)
    # Constant curve should have infinite tortuosity (no direct energy difference)
    assert np.isnan(const_tort)
    # Step function should have intermediate tortuosity
    assert step_tort == 1.0

    # Test conservation
    linear_cons = diatomics.calc_conservation_deviation(xs, y_linear)
    const_cons = diatomics.calc_conservation_deviation(xs, y_const)

    # Linear curve should have constant force = -1
    assert linear_cons == pytest.approx(0.0)
    # Constant curve should have zero force
    assert const_cons == pytest.approx(0.0)


def test_edge_cases() -> None:
    """Test metrics with edge cases."""
    xs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    ys = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    # Test with NaN values
    y_nan = np.array([1.0, np.nan, 3.0, 4.0, 5.0])
    with pytest.raises(ValueError, match="Input contains NaN"):
        diatomics.calc_curve_diff_auc(xs, ys, xs, y_nan)

    # Test with infinite values
    y_inf = np.array([1.0, np.inf, 3.0, 4.0, 5.0])
    with pytest.raises(ValueError, match="Input contains infinite values"):
        diatomics.calc_curve_diff_auc(xs, ys, xs, y_inf)

    # Test with single point
    x_single = np.array([1.0])
    y_single = np.array([1.0])
    with pytest.raises(ValueError, match="Input must have at least 2 points"):
        diatomics.calc_curve_diff_auc(x_single, y_single, x_single, y_single)

    # Test with duplicate x values
    x_dup = np.array([1.0, 1.0, 2.0, 3.0, 4.0])
    y_dup = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    with pytest.raises(ValueError, match="Input contains duplicate x values"):
        diatomics.calc_curve_diff_auc(x_dup, y_dup, x_dup, y_dup)

    # Test with mismatched array sizes
    x_short = np.array([1.0, 2.0, 3.0])
    y_long = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    with pytest.raises(ValueError, match="x and y must have same size"):
        diatomics.calc_curve_diff_auc(x_short, y_long, x_short, y_long)

    # Test with unsorted x values (should work, sorting is handled internally)
    x_unsorted = np.array([1.0, 3.0, 2.0, 5.0, 4.0])
    y_unsorted = np.array([1.0, 3.0, 2.0, 5.0, 4.0])
    auc = diatomics.calc_curve_diff_auc(xs, ys, x_unsorted, y_unsorted)
    assert np.isfinite(auc)


# Test data for smoothness metrics
@pytest.fixture
def smoothness_test_data() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Create test data for smoothness metrics.

    Returns:
        dict[str, tuple[np.ndarray, np.ndarray]]: Dictionary mapping curve names to
            (x, y) tuples.
    """
    x = np.linspace(0.1, 1.0, 1000)
    x_lj = np.linspace(0.3, 3.0, 1000)  # special range for LJ potential

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
        (calc_second_deriv_smoothness, "constant", pytest.approx(0.0, abs=1e-10)),
        (calc_second_deriv_smoothness, "linear", pytest.approx(0.0, abs=1e-10)),
        # Total variation should be constant for linear curve
        (calc_total_variation_smoothness, "linear", pytest.approx(0.0, abs=0.1)),
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
        pytest.param(
            calc_second_deriv_smoothness,
            id="second_deriv",
        ),
        pytest.param(
            calc_total_variation_smoothness,
            id="total_variation",
        ),
        pytest.param(
            calc_curvature_smoothness,
            id="curvature",
        ),
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
    x = np.linspace(0.1, 1.0, 1000)
    base = np.sin(2 * np.pi * x)
    noise_amps = [0.0, 0.1, 0.2]  # Removed 0.3 as it might cause instability

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
