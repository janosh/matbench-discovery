import numpy as np
import pytest

from matbench_discovery import phonons


@pytest.mark.parametrize(
    "freqs, expected",
    [
        (np.array([[0.1, 0.2, 0.3, 0.4]]), False),  # all positive
        (np.array([[0.1, 0.2, -0.3, 0.4]]), True),  # one negative after acoustic
        (np.array([[0.1, 0.2, 0.3], [-0.1, 0.2, 0.3]]), True),  # negative in non-gamma
        (np.array([[-1e-2, -0.1, 0.2, 0.3]]), True),  # below threshold
        (np.full((2, 4), np.nan), True),  # all NaN
    ],
)
def test_check_imaginary_freqs(freqs: np.ndarray, expected: bool) -> None:
    """Test checking for imaginary frequencies."""
    assert phonons.check_imaginary_freqs(freqs) == expected
