"""Tests for geometry optimization metrics (RMSD and symmetry changes vs DFT)."""

from unittest.mock import mock_open, patch

import numpy as np
import pandas as pd
import pytest
from pymatviz.enums import Key

from matbench_discovery.data import make_file_ref
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics import geo_opt

_GEO = "models/alignn/alignn/2026-07-01-geo-opt-symprec=1e-2-moyo=0.4.2.csv.gz"
_GEO_NAN = "models/alignn/alignn/2026-07-02-geo-opt-symprec=1e-2-moyo=0.4.2.csv.gz"


@pytest.mark.parametrize(
    (
        "spg_diffs",
        "n_sym_ops_diffs",
        "expected_decrease",
        "expected_match",
        "expected_increase",
    ),
    [
        # All matches
        ([0, 0, 0], [0, 0, 0], 0.0, 1.0, 0.0),
        # All decreases
        ([-1, -2, -1], [-2, -4, -2], 1.0, 0.0, 0.0),
        # All increases
        ([1, 2, 1], [2, 4, 2], 0.0, 0.0, 1.0),
        # Mixed cases
        ([0, -1, 1], [0, -2, 2], 1 / 3, 1 / 3, 1 / 3),
        # Edge case with zeros in n_sym_ops but non-zero spg
        ([1, -1, 0], [0, 0, 0], 0, 1 / 3, 0),
        # Include some NaN values
        ([0, np.nan, 1], [0, np.nan, 2], 0.0, 0.5, 0.5),
        # NaN in the middle: 5 valid rows with 1 match, 1 decrease, 3 increases
        ([0, 1, np.nan, -1, 2, 3], [0, 2, np.nan, -2, 4, 6], 1 / 5, 1 / 5, 3 / 5),
    ],
)
def test_calc_geo_opt_metrics_parametrized(
    spg_diffs: list[float],
    n_sym_ops_diffs: list[float],
    expected_decrease: float,
    expected_match: float,
    expected_increase: float,
) -> None:
    """Test calc_geo_opt_metrics with various symmetry difference patterns."""
    # RMSD NaN pattern deliberately differs from the spg NaN pattern: n_structures
    # must count valid symmetry rows only and NaN RMSDs are filled with 1.0 (the
    # StructureMatcher stol) before averaging
    rmsds = [0.1] * (len(spg_diffs) - 1) + [np.nan]
    df_geo_opt = pd.DataFrame(
        {
            MbdKey.structure_rmsd_vs_dft: rmsds,
            MbdKey.spg_num_diff: spg_diffs,
            MbdKey.n_sym_ops_diff: n_sym_ops_diffs,
        }
    )

    results = geo_opt.calc_geo_opt_metrics(df_geo_opt)

    assert results[str(Key.symmetry_decrease)] == pytest.approx(expected_decrease)
    assert results[str(Key.symmetry_match)] == pytest.approx(expected_match)
    assert results[str(Key.symmetry_increase)] == pytest.approx(expected_increase)
    # n_structures should be the number of non-NaN spg_diff values
    assert results[str(Key.n_structures)] == np.count_nonzero(pd.notna(spg_diffs))
    expected_rmsd = (0.1 * (len(spg_diffs) - 1) + 1.0) / len(spg_diffs)
    assert results[str(MbdKey.structure_rmsd_vs_dft)] == pytest.approx(expected_rmsd)


@pytest.mark.parametrize(
    ("metrics_data", "expected_block", "analysis_file_path"),
    [
        (
            {
                MbdKey.structure_rmsd_vs_dft: 0.1,
                Key.n_sym_ops_mae: 0.2,
                Key.symmetry_decrease: 0.3,
                Key.symmetry_match: 0.4,
                Key.symmetry_increase: 0.0,
                Key.n_structures: 0,
            },
            {
                Key.rmsd: 0.1,
                Key.n_sym_ops_mae: 0.2,
                Key.symmetry_decrease: 0.3,
                Key.symmetry_match: 0.4,
                Key.symmetry_increase: 0.0,
                Key.n_structures: 0,
                "analysis_file": make_file_ref(_GEO),
            },
            _GEO,
        ),
        (
            {
                MbdKey.structure_rmsd_vs_dft: float("nan"),
                Key.n_sym_ops_mae: float("nan"),
                Key.symmetry_decrease: 0.0,
                Key.symmetry_match: 0.0,
                Key.symmetry_increase: 0.0,
                Key.n_structures: 0,
            },
            {
                Key.rmsd: float("nan"),
                Key.n_sym_ops_mae: float("nan"),
                Key.symmetry_decrease: 0.0,
                Key.symmetry_match: 0.0,
                Key.symmetry_increase: 0.0,
                Key.n_structures: 0,
                "analysis_file": make_file_ref(_GEO_NAN),
            },
            _GEO_NAN,
        ),
    ],
)
def test_write_geo_opt_metrics_to_yaml(
    metrics_data: dict[MbdKey | Key, float],
    expected_block: dict[MbdKey | Key | str, float | str | None],
    analysis_file_path: str,
) -> None:
    """Test saving geometry optimization metrics to YAML files with edge cases."""
    symprec = 1e-2
    symprec_key = "symprec=1e-2"

    # Mock file and YAML operations
    with (
        patch("builtins.open", mock_open()) as mock_file,
        patch("matbench_discovery.data.round_trip_yaml") as mock_yaml,
    ):
        # Configure mock YAML load to return empty dict
        mock_yaml.load.return_value = {}

        # Call the function
        geo_opt.write_metrics_to_yaml(
            pd.DataFrame([metrics_data]), Model.alignn, symprec, analysis_file_path
        )

        # Verify YAML dump was called with expected content
        actual_yaml = mock_yaml.dump.call_args[0][0]

        # Compare metrics while handling NaN values
        actual_block = actual_yaml["metrics"]["geo_opt"][symprec_key]
        assert set(actual_block) == set(expected_block)
        for key, expected_val in expected_block.items():
            value = actual_block[key]
            if key in {"analysis_file", "pred_file", "force_file", "run_info_file"}:
                assert value == expected_val
            else:
                assert value == pytest.approx(expected_val, nan_ok=True)

        # Verify file operations
        mock_file.assert_called()
        mock_yaml.dump.assert_called_once()
