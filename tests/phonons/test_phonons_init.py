"""Tests for phonon helper functions."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pymatviz.enums import Key

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


@pytest.mark.parametrize("id_col", ["mp_id", Key.mat_id])
def test_read_kappa_json_sets_material_id_index(tmp_path: Path, id_col: str) -> None:
    """Kappa JSON files are indexed by material ID, including legacy mp_id."""
    json_path = tmp_path / "kappa.json"
    pd.DataFrame({id_col: ["mp-1", "mp-2"], "kappa": [1.0, 2.0]}).to_json(json_path)

    df_kappa = phonons.read_kappa_json(str(json_path))

    assert list(df_kappa.index) == ["mp-1", "mp-2"]
    assert df_kappa.index.name == Key.mat_id


def test_read_kappa_json_rejects_missing_id_column(tmp_path: Path) -> None:
    """Kappa JSON files without material IDs fail with a clear error."""
    json_path = tmp_path / "kappa.json"
    pd.DataFrame({"kappa": [1.0]}).to_json(json_path)

    with pytest.raises(ValueError, match=r"mp_id.*material_id.*columns=\['kappa'\]"):
        phonons.read_kappa_json(str(json_path))
