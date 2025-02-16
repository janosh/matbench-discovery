import re
from typing import Any

import numpy as np
import pytest

from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics import diatomics
from matbench_discovery.metrics.diatomics import DiatomicCurve

# (ref_curves: {element: (seps, energies)}, pred_curves: {element: (seps, energies)})
PredRefCurves = tuple[dict[str, DiatomicCurve], dict[str, DiatomicCurve]]


@pytest.fixture
def test_data() -> PredRefCurves:
    """Create test data for diatomic curves.

    Returns:
        TestCurves: Reference and predicted curves for each element.
    """
    # Simple test case: H2 molecule with slightly different curves
    xs = np.linspace(0.5, 5.0, 100)
    # Reference: Simple Morse potential
    y_ref = 5 * (1 - np.exp(-2 * (xs - 1.5))) ** 2 - 5
    # Prediction: Slightly perturbed Morse potential
    y_pred = 5.2 * (1 - np.exp(-2.1 * (xs - 1.48))) ** 2 - 5.1

    ref_curves = {"H": (xs, y_ref)}
    pred_curves = {"H": (xs, y_pred)}

    return ref_curves, pred_curves


def test_force_flips(test_data: PredRefCurves) -> None:
    """Test force flips calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    flips = diatomics.calc_force_flips(x_pred, y_pred)
    assert flips == 2

    # Test with different force threshold
    flips_strict = diatomics.calc_force_flips(x_pred, y_pred, threshold=1e-3)
    assert flips_strict >= flips  # Stricter threshold should find more flips


def test_force_total_variation(test_data: PredRefCurves) -> None:
    """Test force total variation calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    total_var = diatomics.calc_force_total_variation(x_pred, y_pred)
    assert isinstance(total_var, float)
    assert total_var >= 0  # Total variation should be non-negative


def test_force_jump(test_data: PredRefCurves) -> None:
    """Test force jump calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    f_jump = diatomics.calc_force_jump(x_pred, y_pred)
    assert isinstance(f_jump, float)
    assert f_jump >= 0  # Force jump should be non-negative


def test_diatomic_curve_metrics(test_data: PredRefCurves) -> None:
    """Test full metrics calculation pipeline."""
    ref_curves, pred_curves = test_data

    # Test with default parameters (no force curves)
    metrics = diatomics.calc_diatomic_curve_metrics(ref_curves, pred_curves)
    assert isinstance(metrics, dict)
    assert "H" in metrics
    metric_keys = [*metrics["H"]]
    assert set(metric_keys) >= {
        MbdKey.norm_auc,
        MbdKey.smoothness,
        MbdKey.tortuosity,
        MbdKey.conservation,
        MbdKey.energy_jump,
        MbdKey.energy_diff_flips,
        MbdKey.energy_grad_norm_max,
    }

    # Test with force curves
    metrics_with_forces = diatomics.calc_diatomic_curve_metrics(
        ref_curves, pred_curves, pred_force_curves=pred_curves
    )
    force_metric_keys = set(metrics_with_forces["H"])
    assert force_metric_keys >= {
        MbdKey.force_total_variation,
        MbdKey.force_jump,
        MbdKey.force_flips,
    }

    # Test with custom parameters
    custom_metrics: dict[str, dict[str, Any]] = {
        MbdKey.norm_auc: {"seps_range": (1.0, 4.0)},
    }
    custom_results = diatomics.calc_diatomic_curve_metrics(
        ref_curves, pred_curves, metrics=custom_metrics
    )
    assert isinstance(custom_results, dict)
    assert "H" in custom_results
    for key in custom_metrics:
        assert key in custom_results["H"]
        assert custom_results["H"][key] != metrics["H"][key], f"{key=}"

    # Test with invalid metric name
    with pytest.raises(
        ValueError, match=re.escape("unknown_metrics={'invalid'}. Valid metrics=")
    ):
        diatomics.calc_diatomic_curve_metrics(
            ref_curves, pred_curves, metrics={"invalid": {}}
        )
