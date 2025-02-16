import numpy as np
import pytest

from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics import diatomics


def test_identical_curves() -> None:
    """Test metrics when comparing a curve to itself."""
    xs = np.linspace(0.5, 5.0, 100)
    ys = 5 * (1 - np.exp(-2 * (xs - 1.5))) ** 2 - 5

    ref_curves = {"H": (xs, ys)}
    pred_curves = {"H": (xs, ys)}

    metrics = diatomics.calc_diatomic_curve_metrics(ref_curves, pred_curves)

    # When curves are identical, certain metrics should be zero or one
    assert metrics["H"][MbdKey.norm_auc] == pytest.approx(0.0)
    assert metrics["H"][MbdKey.energy_mae_vs_ref] == pytest.approx(0.0)
    assert metrics["H"][MbdKey.tortuosity] == pytest.approx(1.0)
    assert metrics["H"][MbdKey.energy_jump] == pytest.approx(0.083043538)


def test_shifted_curves() -> None:
    """Test metrics with shifted curves."""
    xs = np.linspace(0.5, 5.0, 100)
    ys = 5 * (1 - np.exp(-2 * (xs - 1.5))) ** 2 - 5

    # Create shifted curves
    y_shift_up = ys + 1.0  # Constant energy shift
    y_shift_right = 5 * (1 - np.exp(-2 * (xs - 2.0))) ** 2 - 5  # Shifted along x

    ref_curves = {"H": (xs, ys)}
    pred_curves_up = {"H": (xs, y_shift_up)}
    pred_curves_right = {"H": (xs, y_shift_right)}

    # Test energy shift
    metrics_up = diatomics.calc_diatomic_curve_metrics(ref_curves, pred_curves_up)
    assert metrics_up["H"][MbdKey.norm_auc] > 0.0
    # But smoothness should be identical
    assert metrics_up["H"][MbdKey.energy_mae_vs_ref] == pytest.approx(1.0)
    assert metrics_up["H"]["smoothness"] == pytest.approx(578.772692)

    # Test x shift
    metrics_right = diatomics.calc_diatomic_curve_metrics(ref_curves, pred_curves_right)
    assert metrics_right["H"][MbdKey.norm_auc] > 0.0
    assert metrics_right["H"][MbdKey.energy_mae_vs_ref] > 0.0
    assert metrics_right["H"]["smoothness"] == pytest.approx(4561.131837)
