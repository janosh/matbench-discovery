from unittest.mock import mock_open, patch

import numpy as np
import pandas as pd
import pytest
from pymatviz.enums import Key

from matbench_discovery.data import Model
from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics.geo_opt import (
    calc_geo_opt_metrics,
    write_geo_opt_metrics_to_yaml,
)


@pytest.fixture
def df_geo_opt() -> pd.DataFrame:
    """Create a DataFrame simulating geometry optimization data for a single model."""
    return pd.DataFrame(
        {
            MbdKey.structure_rmsd_vs_dft: [0.1, 0.2, 0.3],
            MbdKey.spg_num_diff: [0, -1, 2],
            MbdKey.n_sym_ops_diff: [0, -2, 4],
        }
    )


def test_calc_geo_opt_metrics(df_geo_opt: pd.DataFrame) -> None:
    """Test calc_geo_opt_metrics with a single model."""
    results = calc_geo_opt_metrics(df_geo_opt)

    # Check symmetry change metrics
    assert results[str(Key.symmetry_decrease)] == pytest.approx(1 / 3)
    assert results[str(Key.symmetry_match)] == pytest.approx(1 / 3)
    assert results[str(Key.symmetry_increase)] == pytest.approx(1 / 3)
    assert results[str(Key.n_structures)] == 3


@pytest.mark.parametrize(
    "spg_diffs, n_sym_ops_diffs, expected_decrease, expected_match, expected_increase",
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
    df_geo_opt = pd.DataFrame(
        {
            MbdKey.structure_rmsd_vs_dft: [0.1] * len(spg_diffs),
            MbdKey.spg_num_diff: spg_diffs,
            MbdKey.n_sym_ops_diff: n_sym_ops_diffs,
        }
    )

    results = calc_geo_opt_metrics(df_geo_opt)

    assert results[str(Key.symmetry_decrease)] == pytest.approx(expected_decrease)
    assert results[str(Key.symmetry_match)] == pytest.approx(expected_match)
    assert results[str(Key.symmetry_increase)] == pytest.approx(expected_increase)
    assert results[str(Key.n_structures)] == len(spg_diffs) - np.isnan(spg_diffs).sum()


@pytest.mark.parametrize(
    ("metrics_data", "expected_yaml", "symprec"),
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
                "metrics": {
                    "geo_opt": {
                        "symprec=1e-2": {
                            Key.rmsd: 0.1,
                            Key.n_sym_ops_mae: 0.2,
                            Key.symmetry_decrease: 0.3,
                            Key.symmetry_match: 0.4,
                            Key.symmetry_increase: 0.0,
                            Key.n_structures: 0,
                        }
                    }
                }
            },
            1e-2,
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
                "metrics": {
                    "geo_opt": {
                        "symprec=1e-2": {
                            Key.rmsd: float("nan"),
                            Key.n_sym_ops_mae: float("nan"),
                            Key.symmetry_decrease: 0.0,
                            Key.symmetry_match: 0.0,
                            Key.symmetry_increase: 0.0,
                            Key.n_structures: 0,
                        }
                    }
                }
            },
            1e-2,
        ),
    ],
)
def test_write_geo_opt_metrics_to_yaml(
    metrics_data: dict[str | MbdKey, float],
    expected_yaml: dict[str, dict[str, dict[str, dict[str | MbdKey, float]]]],
    symprec: float,
) -> None:
    """Test saving geometry optimization metrics to YAML files with edge cases."""
    symprec_key = f"{symprec=:.0e}".replace("e-0", "e-")

    # Mock the Model class and file operations
    with (
        patch("matbench_discovery.metrics.geo_opt.Model") as mock_model,
        patch("builtins.open", mock_open()) as mock_file,
    ):
        # Configure mock model
        mock_model.from_label.return_value.label = "test_model"
        mock_model.from_label.return_value.yaml_path = "mock_path/test_model.yml"

        # Mock the YAML operations
        with patch("matbench_discovery.metrics.geo_opt.round_trip_yaml") as mock_yaml:
            # Configure mock YAML load to return empty dict
            mock_yaml.load.return_value = {}

            # Call the function
            write_geo_opt_metrics_to_yaml(metrics_data, Model.alignn, symprec)

            # Verify YAML dump was called with expected content
            actual_yaml = mock_yaml.dump.call_args[0][0]

            # Compare metrics while handling NaN values
            for key, value in actual_yaml["metrics"]["geo_opt"][symprec_key].items():
                expected_value = expected_yaml["metrics"]["geo_opt"][symprec_key][key]
                if isinstance(value, float) and np.isnan(value):
                    assert np.isnan(expected_value)
                else:
                    assert value == pytest.approx(expected_value)

            # Verify file operations
            mock_file.assert_called()
            mock_yaml.dump.assert_called_once()
