from unittest.mock import mock_open, patch

import numpy as np
import pandas as pd
import pytest
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey
from matbench_discovery.metrics.geo_opt import (
    analyze_symmetry_changes,
    write_geo_opt_metrics_to_yaml,
)

model1, model2, model3 = "model1", "model2", "model3"


@pytest.fixture
def df_sym() -> pd.DataFrame:
    """Create a multi-index DataFrame simulating symmetry data."""

    data = {
        (model1, MbdKey.structure_rmsd_vs_dft): [0.1, 0.2, 0.3],
        (model1, MbdKey.spg_num_diff): [0, -1, 2],
        (model1, MbdKey.n_sym_ops_diff): [0, -2, 4],
        (model2, MbdKey.structure_rmsd_vs_dft): [0.2, 0.3, 0.4],
        (model2, MbdKey.spg_num_diff): [1, 0, -1],
        (model2, MbdKey.n_sym_ops_diff): [2, 0, -2],
        (Key.dft.label, MbdKey.structure_rmsd_vs_dft): [0, 0, 0],
        (Key.dft.label, MbdKey.spg_num_diff): [0, 0, 0],
        (Key.dft.label, MbdKey.n_sym_ops_diff): [0, 0, 0],
    }
    return pd.DataFrame(data)


def test_analyze_symmetry_changes(df_sym: pd.DataFrame) -> None:
    """Test analyze_symmetry_changes with multiple models."""
    results = analyze_symmetry_changes(df_sym)

    # Check results for model1
    sym_change_cols = [Key.symmetry_decrease, Key.symmetry_match, Key.symmetry_increase]
    assert results.loc[model1, sym_change_cols].to_list() == [1 / 3, 1 / 3, 1 / 3]
    assert results.loc[model1, Key.n_structs] == 3

    # Check results for model2
    assert results.loc[model2, sym_change_cols].to_list() == [1 / 3, 1 / 3, 1 / 3]
    assert results.loc[model2, Key.n_structs] == 3

    # Check that DFT is not included in results
    assert Key.dft.label not in results.index


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
def test_analyze_symmetry_changes_parametrized(
    spg_diffs: list[float],
    n_sym_ops_diffs: list[float],
    expected_decrease: float,
    expected_match: float,
    expected_increase: float,
) -> None:
    """Test analyze_symmetry_changes with various symmetry difference patterns."""
    model = "test_model"
    df_sym = pd.DataFrame(
        {
            (model, MbdKey.structure_rmsd_vs_dft): [0.1] * len(spg_diffs),
            (model, MbdKey.spg_num_diff): spg_diffs,
            (model, MbdKey.n_sym_ops_diff): n_sym_ops_diffs,
        }
    )

    results = analyze_symmetry_changes(df_sym)

    assert results.loc[model, Key.symmetry_decrease] == pytest.approx(expected_decrease)
    assert results.loc[model, Key.symmetry_match] == pytest.approx(expected_match)
    assert results.loc[model, Key.symmetry_increase] == pytest.approx(expected_increase)
    assert (
        results.loc[model, Key.n_structs] == len(spg_diffs) - np.isnan(spg_diffs).sum()
    )


@pytest.fixture
def df_sym_changes() -> pd.DataFrame:
    """Create a DataFrame with symmetry change statistics."""
    return pd.DataFrame(
        {
            Key.symmetry_decrease: [0.2, 0.3],
            Key.symmetry_match: [0.5, 0.4],
            Key.symmetry_increase: [0.3, 0.3],
            Key.n_structs: [100, 100],
        },
        index=["model1", "model2"],
    )


sym_changes_cols = [
    Key.symmetry_decrease,
    Key.symmetry_match,
    Key.symmetry_increase,
    Key.n_structs,
]


@pytest.mark.parametrize(
    "df_sym_data, df_sym_changes_data, expected_yaml",
    [
        # Test case 1: Normal case with valid metrics
        (
            {  # df_sym data
                (model1, MbdKey.structure_rmsd_vs_dft): [0.1, 0.2, 0.3],
                (model1, MbdKey.spg_num_diff): [0, 1, -1],
                (model1, MbdKey.n_sym_ops_diff): [0, 2, -2],
            },
            # df_sym_changes data
            model1_sym_changes := dict(zip(sym_changes_cols, [0.33, 0.33, 0.34, 3])),
            # expected yaml content
            {"metrics": {"geo_opt": {**model1_sym_changes, Key.rmsd: 0.2}}},
        ),
        # Test case 2: Edge case with NaN values
        (
            {
                (model2, MbdKey.structure_rmsd_vs_dft): [np.nan, np.nan],
                (model2, MbdKey.spg_num_diff): [np.nan, np.nan],
                (model2, MbdKey.n_sym_ops_diff): [np.nan, np.nan],
            },
            model2_sym_changes := dict(zip(sym_changes_cols, [0.0, 0.0, 0.0, 0])),
            {"metrics": {"geo_opt": {**model2_sym_changes, Key.rmsd: float("nan")}}},
        ),
        # Test case 3: Empty data
        (
            {
                (model3, MbdKey.structure_rmsd_vs_dft): [],
                (model3, MbdKey.spg_num_diff): [],
                (model3, MbdKey.n_sym_ops_diff): [],
            },
            model3_sym_changes := dict(zip(sym_changes_cols, [0.0, 0.0, 0.0, 0])),
            {"metrics": {"geo_opt": {**model3_sym_changes, Key.rmsd: float("nan")}}},
        ),
    ],
)
def test_write_geo_opt_metrics_to_yaml(
    df_sym_data: dict[tuple[str, MbdKey], list[float]],
    df_sym_changes_data: dict[str, list[float]],
    expected_yaml: dict[str, dict[str, dict[str, float | int]]],
) -> None:
    """Test saving geometry optimization metrics to YAML files with edge cases."""
    # Create test DataFrames
    df_sym = pd.DataFrame(df_sym_data)
    # Set the column names for the MultiIndex
    df_sym.columns.names = ["model", MbdKey.sym_prop]

    model_name = df_sym.columns.levels[0][0]
    df_sym_changes = pd.DataFrame([df_sym_changes_data], index=[model_name])

    # Mock the Model class and file operations
    with (
        patch("matbench_discovery.metrics.geo_opt.Model") as mock_model,
        patch("builtins.open", mock_open()) as mock_file,
    ):
        # Configure mock model
        mock_model.from_label.return_value.label = model_name
        mock_model.from_label.return_value.yaml_path = f"mock_path/{model_name}.yml"

        # Mock the YAML operations
        with patch("matbench_discovery.metrics.geo_opt.round_trip_yaml") as mock_yaml:
            # Configure mock YAML load to return empty dict
            mock_yaml.load.return_value = {}

            # Call the function
            write_geo_opt_metrics_to_yaml(df_sym, df_sym_changes)

            # Verify YAML dump was called with expected content
            actual_yaml = mock_yaml.dump.call_args[0][0]

            # Compare metrics while handling NaN values
            for key, value in actual_yaml["metrics"]["geo_opt"].items():
                expected_value = expected_yaml["metrics"]["geo_opt"][key]
                if isinstance(value, float) and np.isnan(value):
                    assert np.isnan(expected_value)
                else:
                    assert value == pytest.approx(expected_value)

            # Verify file operations
            mock_file.assert_called()
            mock_yaml.dump.assert_called_once()
