from collections.abc import Callable

import numpy as np
import pytest

from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics import diatomics


def base_curve(xs: np.ndarray) -> np.ndarray:
    """Generate base test curve."""
    return 5 * (1 - np.exp(-2 * (xs - 1.5))) ** 2 - 5


@pytest.mark.parametrize(
    "curve_func,expected",
    [
        (  # No modification (identical curves ref_ys == pred_ys)
            lambda xs: base_curve(xs),
            (0.0, 0.0, 1.0, 0.083043538, 578.772692),
        ),
        (  # Vertical shift (+1.0)
            lambda xs: base_curve(xs) + 1.0,
            (0.00489955, 1.0, 1.0, 0.083043538, 578.772692),
        ),
        (  # Horizontal shift (+0.5)
            lambda xs: 5 * (1 - np.exp(-2 * (xs - 2.0))) ** 2 - 5,
            (0.411922297, 84.797947, 1.0, 0.083043538, 4561.131837),
        ),
    ],
)
def test_curve_shifts(
    curve_func: Callable[[np.ndarray], np.ndarray], expected: tuple[float, ...]
) -> None:
    """Test metrics for various curve modifications.

    Args:
        curve_func (Callable): Function that takes x values and returns y values
        expected (tuple[float, ...]): Expected values for metrics in order:
            (norm_auc, energy_mae_vs_ref, tortuosity, energy_jump, smoothness)
    """
    xs = np.linspace(0.5, 5.0, 100)
    ref_ys = base_curve(xs)
    pred_ys = curve_func(xs)

    ref_curves = {"H": (xs, ref_ys)}
    pred_curves = {"H": (xs, pred_ys)}

    metrics = diatomics.calc_diatomic_curve_metrics(ref_curves, pred_curves)

    # Check metrics
    metric_keys = (
        MbdKey.norm_auc,
        MbdKey.energy_mae_vs_ref,
        MbdKey.tortuosity,
        MbdKey.energy_jump,
        MbdKey.smoothness,
    )
    for metric_key, expect in zip(metric_keys, expected, strict=True):
        actual = metrics["H"][metric_key]
        assert actual == pytest.approx(expect), (
            f"{metric_key=} expected {expect}, got {actual}"
        )
