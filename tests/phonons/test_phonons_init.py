"""Tests for phonon helper functions."""

import json
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


def test_read_kappa_json_sets_material_id_index(tmp_path: Path) -> None:
    """Kappa JSON is indexed by material ID; only flagged records drop conductivity
    and NaN-filled optional flags are derived as if absent.
    """
    json_path = tmp_path / "kappa.json"
    tensor = np.eye(3).tolist()
    # the shared runner writes conductivity_skipped only on skipped records; mp-1
    # lacks broken_symmetry, which the DataFrame fills with NaN from mp-3's column
    records = [
        {
            Key.mat_id: "mp-1",
            "kappa_tot_rta": [tensor],
            "errors": [],
            Key.init_spg_num: 225,
            Key.final_spg_num: 221,
        },
        {Key.mat_id: "mp-2", "kappa_tot_rta": [tensor], "conductivity_skipped": True},
        {Key.mat_id: "mp-3", "kappa_tot_rta": [tensor], "broken_symmetry": False},
    ]
    json_path.write_text(json.dumps(records))

    df_kappa = phonons.read_kappa_json(str(json_path))

    assert list(df_kappa.index) == ["mp-1", "mp-2", "mp-3"]
    assert df_kappa.index.name == Key.mat_id
    np.testing.assert_array_equal(df_kappa.loc["mp-1", "kappa_tot_rta"], [tensor])
    assert np.isnan(df_kappa.loc["mp-2", "kappa_tot_rta"])
    assert df_kappa["broken_symmetry"].tolist() == [True, False, False]


@pytest.mark.parametrize("ids", [None, [""], [None], [1], ["mp-1", "mp-1"]])
def test_read_kappa_json_rejects_invalid_ids(tmp_path: Path, ids: list | None) -> None:
    """Missing, empty, non-string, duplicate, and legacy-only IDs fail early."""
    json_path = tmp_path / "kappa.json"
    pd.DataFrame({"mp_id": ["mp-1"]} if ids is None else {Key.mat_id: ids}).to_json(
        json_path
    )

    with pytest.raises(ValueError, match=r"material_id|duplicate"):
        phonons.read_kappa_json(str(json_path))
