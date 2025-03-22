from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import yaml
from pymatviz.enums import Key

from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics.geo_opt import (
    calc_geo_opt_metrics,
    write_metrics_to_yaml,
)


@pytest.fixture
def df_geo_opt() -> pd.DataFrame:
    """Create a test DataFrame with geometry optimization metrics."""
    return pd.DataFrame(
        {
            MbdKey.spg_num_diff: [0, -1, 1],
            MbdKey.n_sym_ops_diff: [0, -2, 2],
            MbdKey.structure_rmsd_vs_dft: [0.1, 0.2, 0.3],
        }
    )


@pytest.mark.parametrize(
    "test_data, expected",
    [
        (  # Basic test case with both symmetry and distance
            {
                MbdKey.spg_num_diff: [0, -1, 1],
                MbdKey.n_sym_ops_diff: [0, -2, 2],
                MbdKey.structure_rmsd_vs_dft: [0.1, 0.2, 0.3],
            },
            {
                str(Key.n_structures): 3,
                str(MbdKey.structure_rmsd_vs_dft): 0.2,
                str(Key.n_sym_ops_mae): 4 / 3,
                str(Key.symmetry_match): 1 / 3,
                str(Key.symmetry_decrease): 1 / 3,
                str(Key.symmetry_increase): 1 / 3,
            },
        ),
        (  # Distance only metrics
            # NaN replaced with 1.0
            {MbdKey.structure_rmsd_vs_dft: [0.1, 0.2, np.nan, 0.3]},
            {str(Key.n_structures): 4, str(MbdKey.structure_rmsd_vs_dft): 0.4},
            # Symmetry metrics should be absent
        ),
        (  # Symmetry only metrics
            {
                MbdKey.spg_num_diff: [0, 0, 0],
                MbdKey.n_sym_ops_diff: [0, 0, 0],
            },
            {
                str(Key.n_structures): 3,
                str(Key.symmetry_match): 1.0,
                str(Key.symmetry_decrease): 0.0,
                str(Key.symmetry_increase): 0.0,
                str(Key.n_sym_ops_mae): 0.0,
            },
        ),
        (  # With NaN values
            {
                MbdKey.spg_num_diff: [0, np.nan, 1],
                MbdKey.n_sym_ops_diff: [0, np.nan, 2],
                MbdKey.structure_rmsd_vs_dft: [0.1, 0.2, np.nan],
            },
            {
                str(Key.n_structures): 3,  # this is total row count, not valid ones
                str(MbdKey.structure_rmsd_vs_dft): 0.433,  # (0.1 + 0.2 + 1.0) / 3
                str(Key.symmetry_match): 0.5,  # 1/2
                str(Key.symmetry_decrease): 0.0,
                str(Key.symmetry_increase): 0.5,  # 1/2
            },
        ),
    ],
)
def test_calc_geo_opt_metrics(
    test_data: dict[str, list[float | None]], expected: dict[str, float]
) -> None:
    """Test that geo_opt metrics are correctly calculated for different inputs."""
    test_df = pd.DataFrame(test_data)
    metrics = calc_geo_opt_metrics(test_df)

    # Check all expected metrics exist with correct values
    for key, value in expected.items():
        if isinstance(value, float):
            assert metrics[key] == pytest.approx(value, abs=1e-3)
        else:
            assert metrics[key] == value

    # Check symmetry metrics are only present when expected
    symmetry_keys = [
        str(Key.symmetry_match),
        str(Key.symmetry_decrease),
        str(Key.symmetry_increase),
        str(Key.n_sym_ops_mae),
    ]

    has_symmetry_data = (
        MbdKey.spg_num_diff in test_data and MbdKey.n_sym_ops_diff in test_data
    )

    for key in symmetry_keys:
        if has_symmetry_data:
            assert key in metrics
        else:
            assert key not in metrics


@pytest.mark.parametrize(
    "input_data, section, symprec, file_path, expected_checks, error_expected",
    [
        (  # Test DataFrame input with symmetry section
            "df_geo_opt",  # Use fixture
            "symmetry",
            1e-2,
            "test/path/df_input.csv",
            {
                "metrics.geo_opt.symmetry.symprec=1e-2": {
                    # (value, tolerance) for approx comparison
                    str(Key.symmetry_match): (1 / 3, 1e-4),
                    str(Key.symmetry_decrease): (1 / 3, 1e-4),
                    str(Key.symmetry_increase): (1 / 3, 1e-4),
                    str(Key.n_sym_ops_mae): (4 / 3, 1e-4),
                },
                "metrics.geo_opt.analysis_file": "test/path/df_input.csv",
            },
            False,
        ),
        (  # Test dict input with symmetry section
            {
                str(MbdKey.structure_rmsd_vs_dft): 0.1,
                str(Key.n_sym_ops_mae): 0.2,
                str(Key.symmetry_decrease): 0.3,
                str(Key.symmetry_match): 0.4,
                str(Key.symmetry_increase): 0.5,
                str(Key.n_structures): 10,
            },
            "symmetry",
            1e-2,
            "test/path/symmetry.csv",
            {
                "metrics.geo_opt.symmetry.symprec=1e-2": {
                    str(Key.n_sym_ops_mae): 0.2,
                    str(Key.symmetry_decrease): 0.3,
                    str(Key.symmetry_match): 0.4,
                    str(Key.symmetry_increase): 0.5,
                    str(Key.rmsd): None,  # Should not be present
                },
                "metrics.geo_opt.analysis_file": "test/path/symmetry.csv",
            },
            False,
        ),
        (  # Test dict input with distance section
            {
                str(MbdKey.structure_rmsd_vs_dft): 0.15,
                str(Key.n_structures): 20,
            },
            "distance",
            None,  # symprec not needed for distance
            "test/path/distance.csv",
            {
                "metrics.geo_opt.distance": {
                    str(Key.rmsd): 0.15,
                    str(Key.n_sym_ops_mae): None,  # Should not be present
                },
                "metrics.geo_opt.n_structures": 20,
            },
            False,
        ),
        (  # Test error case: symmetry section without symprec
            {
                str(Key.n_sym_ops_mae): 0.2,
                str(Key.symmetry_decrease): 0.3,
            },
            "symmetry",
            None,  # Missing symprec should cause error
            "test/path/error.csv",
            {},
            True,
        ),
    ],
)
def test_write_metrics_to_yaml(
    input_data: dict[str, float] | str,
    section: str,
    symprec: float | None,
    file_path: str,
    expected_checks: dict[str, Any],
    error_expected: bool,
    tmp_path: Path,
    df_geo_opt: pd.DataFrame,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test that metrics are correctly written to YAML under different scenarios."""
    yaml_path = f"{tmp_path}/test_model.yml"

    # Create an empty YAML file
    with open(yaml_path, mode="w") as file:
        file.write("metrics: {}\n")

    # Handle fixture reference
    if input_data == "df_geo_opt":
        input_data = df_geo_opt

    # Use a real Model enum (alchembert) with a patched yaml_path
    test_model = Model.alchembert

    def mock_yaml_path(_self: Model) -> str:
        """Mock implementation that returns the test yaml path."""
        return yaml_path

    # Apply the patch using monkeypatch
    monkeypatch.setattr(Model, "yaml_path", property(mock_yaml_path))

    # Test function call
    if error_expected:
        with pytest.raises(ValueError, match="symprec must be provided"):
            write_metrics_to_yaml(input_data, test_model, file_path, section, symprec)
        return

    # Normal execution
    write_metrics_to_yaml(input_data, test_model, file_path, section, symprec)

    # Load and check the file
    with open(yaml_path) as file:
        data = yaml.safe_load(file)

    # Verify expected values
    for path, checks in expected_checks.items():
        # Navigate to the specified path in the YAML
        path_parts = path.split(".")
        current = data
        for part in path_parts:
            if part in current:
                current = current[part]
            else:
                pytest.fail(f"Path {path} not found in YAML data")

        # Handle different types of checks
        if not isinstance(checks, dict):
            # direct value comparison
            # (e.g. "metrics.geo_opt.analysis_file": "test/path.csv")
            assert current == checks
            continue

        # Check values at this path for dictionary checks
        for key, expected_value in checks.items():
            if expected_value is None:
                # Check that key is not present
                assert key not in current
            elif isinstance(expected_value, tuple) and len(expected_value) == 2:
                # Handle approx comparison with (value, tolerance)
                value, tolerance = expected_value
                assert float(current[key]) == pytest.approx(value, abs=tolerance)
            else:
                assert current[key] == expected_value
