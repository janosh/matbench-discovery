import re
from typing import Any

import numpy as np
import pytest

from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics import diatomics


def test_force_flips(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test force flips calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    flips = diatomics.calc_force_flips(x_pred, y_pred)
    assert flips == 2

    # Test with different force threshold
    flips_strict = diatomics.calc_force_flips(x_pred, y_pred, threshold=1e-3)
    assert flips_strict >= flips  # Stricter threshold should find more flips


def test_force_total_variation(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test force total variation calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    total_var = diatomics.calc_force_total_variation(x_pred, y_pred)
    assert isinstance(total_var, float)
    assert total_var >= 0  # Total variation should be non-negative


def test_force_jump(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test force jump calculation."""
    _, pred_curves = test_data
    x_pred, y_pred = pred_curves["H"]

    # Test with default parameters
    fjump = diatomics.calc_force_jump(x_pred, y_pred)
    assert isinstance(fjump, float)
    assert fjump >= 0  # Force jump should be non-negative


def test_diatomic_curve_metrics(
    test_data: tuple[
        dict[str, tuple[np.ndarray, np.ndarray]],
        dict[str, tuple[np.ndarray, np.ndarray]],
    ],
) -> None:
    """Test full metrics calculation pipeline."""
    ref_curves, pred_curves = test_data

    # Test with default parameters (no force curves)
    metrics = diatomics.calc_diatomic_curve_metrics(ref_curves, pred_curves)
    assert isinstance(metrics, dict)
    assert "H" in metrics
    metric_keys = [*metrics["H"]]
    for metric in [
        MbdKey.norm_auc,
        MbdKey.smoothness,
        MbdKey.tortuosity,
        MbdKey.conservation,
        MbdKey.energy_jump,
        MbdKey.energy_diff_flips,
        MbdKey.energy_grad_norm_max,
    ]:
        assert metric in metric_keys, f"{metric=} not in {metric_keys=}"

    # Test with force curves
    metrics_with_forces = diatomics.calc_diatomic_curve_metrics(
        ref_curves, pred_curves, pred_force_curves=pred_curves
    )
    force_metric_keys = [*metrics_with_forces["H"]]
    for metric in [
        MbdKey.force_total_variation,
        MbdKey.force_jump,
    ]:
        assert metric in force_metric_keys, f"{metric=} not in {force_metric_keys=}"

    # Test with custom parameters
    custom_metrics: dict[str, dict[str, Any]] = {
        MbdKey.norm_auc: {"seps_range": (1.0, 4.0)},
    }
    custom_results = diatomics.calc_diatomic_curve_metrics(
        ref_curves,
        pred_curves,
        metrics=custom_metrics,
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
