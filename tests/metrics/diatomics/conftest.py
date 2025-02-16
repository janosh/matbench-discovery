import numpy as np
import pytest

from matbench_discovery.metrics.diatomics import DiatomicCurve


@pytest.fixture
def test_data() -> tuple[dict[str, DiatomicCurve], dict[str, DiatomicCurve]]:
    """Create test data for diatomic curves.

    Returns:
        tuple[dict[str, DiatomicCurve], dict[str, DiatomicCurve]]: Reference and
            predicted curves for each element.
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
