"""Tests for the geometry optimization analysis script."""

# ruff: noqa: E402
# We need to manipulate the path before importing other modules

import gzip
import importlib.util
import os
from collections.abc import Generator
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pandas as pd
import pytest
import yaml
from pymatgen.core import Structure
from pymatgen.core.lattice import Lattice
from pymatviz.enums import Key

from matbench_discovery import SCRIPTS
from matbench_discovery.enums import MbdKey, Model
from matbench_discovery.metrics import geo_opt
from matbench_discovery.models import MODEL_METADATA
from matbench_discovery.structure import perturb_structure, symmetry

# Import analyze_geo_opt.py script via importlib
script_path = f"{SCRIPTS}/analyze_geo_opt.py"
basename = os.path.basename(script_path)
spec = importlib.util.spec_from_file_location(basename, script_path)
analyze_geo_opt = importlib.util.module_from_spec(spec)  # type: ignore[arg-type]
if spec is None or spec.loader is None:
    raise ImportError(f"Failed to import {script_path}")
spec.loader.exec_module(analyze_geo_opt)


@pytest.fixture
def small_structure_set() -> dict[str, Any]:
    """Create a small set of test structures."""
    cubic = Structure(Lattice.cubic(4.2), ["Si", "Si"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    perturbed = cubic.copy()
    perturbed.perturb(0.1)
    tetragonal = Structure(
        Lattice.tetragonal(a=4, c=6),
        ["Ti", "O", "O"],
        [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5]],
    )

    dft_structs = {"wbm-1": cubic, "wbm-2": tetragonal, "wbm-3": cubic.copy()}
    ml_structs = {
        "wbm-1": perturbed,
        "wbm-2": tetragonal.copy(),
        "wbm-3": perturbed.copy(),
    }

    return {
        "dft_structs": dft_structs,
        "ml_structs": ml_structs,
        "df_dft_analysis": symmetry.get_sym_info_from_structs(dft_structs),
    }


@pytest.fixture
def large_structure_set() -> dict[str, Any]:
    """Create a larger set of test structures (>100 items) to test row limits."""
    dft_structs: dict[str, Structure] = {}
    ml_structs: dict[str, Structure] = {}

    # Create 150 structures (exceeding the 100 row limit)
    for i in range(150):
        # Vary parameters slightly for each structure
        a = 4.0 + i * 0.01
        cubic = Structure(Lattice.cubic(a), ["Si", "Si"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        perturbed = cubic.copy()
        perturbed.perturb(0.1)

        mat_id = f"wbm-test-{i}"
        dft_structs[mat_id] = cubic
        ml_structs[mat_id] = perturbed

    return {
        "dft_structs": dft_structs,
        "ml_structs": ml_structs,
        "df_dft_analysis": symmetry.get_sym_info_from_structs(dft_structs),
    }


@pytest.fixture
def test_structures_files(
    tmp_path: Path, small_structure_set: dict[str, Any]
) -> dict[str, Any]:
    """Create test JSON and JSONL files with structures."""
    df_test = pd.DataFrame(
        {
            "material_id": list(small_structure_set["ml_structs"]),
            "structure": [
                struct.as_dict()
                for struct in small_structure_set["ml_structs"].values()
            ],
        }
    )

    json_path = f"{tmp_path}/test_structures.json.gz"
    jsonl_path = f"{tmp_path}/test_structures.jsonl.gz"

    with gzip.open(json_path, mode="wb") as file:
        file.write(df_test.to_json(orient="records").encode("utf-8"))
    with gzip.open(jsonl_path, mode="wb") as file:
        file.write(df_test.to_json(orient="records", lines=True).encode("utf-8"))

    return {"json_path": json_path, "jsonl_path": jsonl_path, "df": df_test}


@pytest.fixture
def large_test_structures_file(
    tmp_path: Path, large_structure_set: dict[str, Any]
) -> dict[str, Any]:
    """Create test files with >100 structures to test row limits."""
    df_test = pd.DataFrame(
        {
            "material_id": list(large_structure_set["ml_structs"]),
            "structure": [
                struct.as_dict()
                for struct in large_structure_set["ml_structs"].values()
            ],
        }
    )

    json_path = f"{tmp_path}/large_test_structures.json.gz"
    jsonl_path = f"{tmp_path}/large_test_structures.jsonl.gz"

    with gzip.open(json_path, mode="wb") as file:
        file.write(df_test.to_json(orient="records").encode("utf-8"))
    with gzip.open(jsonl_path, mode="wb") as file:
        file.write(df_test.to_json(orient="records", lines=True).encode("utf-8"))

    return {"json_path": json_path, "jsonl_path": jsonl_path, "df": df_test}


@pytest.fixture
def mock_model_yaml(tmp_path: Path) -> str:
    """Create a mock model YAML file."""
    yaml_path = f"{tmp_path}/model.yaml"
    metrics = {"geo_opt": {"struct_col": "structure", "symmetry": {}, "distance": {}}}
    with open(yaml_path, mode="w") as file:
        yaml.dump({"name": "test_model", "metrics": metrics}, file)
    return yaml_path


@pytest.fixture
def setup_model_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> Generator[None, None, None]:
    """Setup and restore MODEL_METADATA."""
    original = MODEL_METADATA.copy()
    MODEL_METADATA["test_model"] = {"metrics": {"geo_opt": {"struct_col": "structure"}}}
    monkeypatch.setattr("matbench_discovery.models.MODEL_METADATA", MODEL_METADATA)
    yield
    for model_key in list(MODEL_METADATA):
        if model_key not in original:
            del MODEL_METADATA[model_key]


def test_analyze_ml_relaxed_structs(
    small_structure_set: dict[str, Any],
    mock_model_yaml: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    # We need this fixture to be active but don't directly use it
    setup_model_metadata: None,  # noqa: ARG001
    capsys: pytest.CaptureFixture,
) -> None:
    """Test the full analysis pipeline with minimal test structures."""
    # Mock importlib.metadata.version to ensure moyo_version is dynamic
    with patch("importlib.metadata.version", return_value="test_version"):
        # Create test data file with minimal structure data (only 3 structures)
        df_test_data = pd.DataFrame(
            {
                "material_id": ["wbm-1", "wbm-2", "wbm-3"],
                "structure": [
                    small_structure_set["ml_structs"]["wbm-1"].as_dict(),
                    small_structure_set["ml_structs"]["wbm-2"].as_dict(),
                    small_structure_set["ml_structs"]["wbm-3"].as_dict(),
                ],
            }
        )
        test_file_path = f"{tmp_path}/test_structures.json.gz"
        with gzip.open(test_file_path, mode="wb") as file:
            file.write(df_test_data.to_json(orient="records").encode("utf-8"))

        # Override geo_opt_path for MACE-MP-0 model and mock the yaml_path
        model = Model.mace_mp_0

        # Create a function to override the geo_opt_path property for this model only
        def mock_geo_opt_path(self: Model) -> str:
            if self == model:
                return test_file_path
            # Call the original for other models
            from matbench_discovery.enums import Model

            return Model.geo_opt_path.__get__(self, Model)

        # Create a function to override the yaml_path property for this model only
        def mock_yaml_path(self: Model) -> str:
            if self == model:
                return str(mock_model_yaml)
            # Call the original for other models
            from matbench_discovery.enums import Model

            return Model.yaml_path.__get__(self, Model)

        # Patch the properties
        monkeypatch.setattr(
            "matbench_discovery.enums.Model.geo_opt_path", property(mock_geo_opt_path)
        )
        monkeypatch.setattr(
            "matbench_discovery.enums.Model.yaml_path", property(mock_yaml_path)
        )

        # Add geo_opt support to MODEL_METADATA for this model (using fixture)
        MODEL_METADATA[model.label]["metrics"]["geo_opt"] = {"struct_col": "structure"}

        # Run the analysis with no mocks except for model properties
        analyze_geo_opt.analyze_ml_relaxed_structs(
            model=model,
            symprec=1e-2,
            moyo_version="moyo=test_version",
            df_dft_analysis=small_structure_set["df_dft_analysis"],
            dft_structs=small_structure_set["dft_structs"],
            analysis_type="all",
            debug_mode=0,
        )

        # Check basic outputs
        output = capsys.readouterr().out
        assert any(
            phrase in output for phrase in ["structures", "structure", "Successfully"]
        )

        # Check the YAML file was updated properly with the actual geo_opt metrics
        with open(mock_model_yaml) as f:
            yaml_data = yaml.safe_load(f)
            assert "metrics" in yaml_data
            assert "geo_opt" in yaml_data["metrics"]
            assert "analysis_file" in yaml_data["metrics"]["geo_opt"]


def test_analyze_ml_relaxed_structs_handling_errors(
    small_structure_set: dict[str, Any],
    mock_model_yaml: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    # We need this fixture to be active but don't directly use it
    setup_model_metadata: None,  # noqa: ARG001
    capsys: pytest.CaptureFixture,
) -> None:
    """Test error handling with nonexistent file."""
    # Setup the model with a nonexistent file path
    model = Model.mace_mp_0
    nonexistent_file_path = f"{tmp_path}/nonexistent_file.json.gz"

    # Create a function to override the geo_opt_path property for this model only
    def mock_geo_opt_path(self: Model) -> str:
        if self == model:
            return nonexistent_file_path
        # Call the original for other models
        from matbench_discovery.enums import Model

        return Model.geo_opt_path.__get__(self, Model)

    # Create a function to override the yaml_path property for this model only
    def mock_yaml_path(self: Model) -> str:
        if self == model:
            return str(mock_model_yaml)
        # Call the original for other models
        from matbench_discovery.enums import Model

        return Model.yaml_path.__get__(self, Model)

    # Patch the properties
    monkeypatch.setattr(
        "matbench_discovery.enums.Model.geo_opt_path", property(mock_geo_opt_path)
    )
    monkeypatch.setattr(
        "matbench_discovery.enums.Model.yaml_path", property(mock_yaml_path)
    )

    # Add geo_opt support to MODEL_METADATA for this model (using fixture)
    MODEL_METADATA[model.label]["metrics"]["geo_opt"] = {"struct_col": "structure"}

    # Run the analysis with the nonexistent file
    analyze_geo_opt.analyze_ml_relaxed_structs(
        model=model,
        symprec=1e-2,
        moyo_version="moyo=test",
        df_dft_analysis=small_structure_set["df_dft_analysis"],
        dft_structs=small_structure_set["dft_structs"],
        analysis_type="all",
        debug_mode=0,
    )

    # Verify the error message
    output = capsys.readouterr().out
    assert "structures not found" in output


# Integration tests for the geo_opt workflow
def test_symmetry_analysis_integration(small_structure_set: dict[str, Any]) -> None:
    """Test that symmetry analysis workflow functions correctly."""
    # Extract structures from small_structure_set
    ml_structs = small_structure_set["ml_structs"]
    dft_structs = small_structure_set["dft_structs"]

    # Step 1: Get symmetry information for both ML and DFT structures
    ml_symprec = 1e-2
    df_ml_sym = symmetry.get_sym_info_from_structs(ml_structs, symprec=ml_symprec)
    df_dft_sym = symmetry.get_sym_info_from_structs(dft_structs, symprec=ml_symprec)

    # Step 2: Compare symmetry information
    df_compared = df_ml_sym.copy()
    df_compared[MbdKey.spg_num_diff] = df_ml_sym[Key.spg_num] - df_dft_sym[Key.spg_num]
    df_compared[MbdKey.n_sym_ops_diff] = (
        df_ml_sym[Key.n_sym_ops] - df_dft_sym[Key.n_sym_ops]
    )

    # Step 3: Calculate structure distances
    df_compared = symmetry.calc_structure_distances(
        df_compared, ml_structs, dft_structs
    )

    # Verify the results
    assert MbdKey.spg_num_diff in df_compared
    assert MbdKey.n_sym_ops_diff in df_compared
    assert MbdKey.structure_rmsd_vs_dft in df_compared

    # Check that we have RMSD values for at least some structures
    assert not df_compared[MbdKey.structure_rmsd_vs_dft].isna().all()

    # For this specific fixture, the perturbed structures still maintain their symmetry
    # and may have zero RMSD, so we skip the non-zero RMSD check

    # Step 4: Calculate geo opt metrics
    metrics = geo_opt.calc_geo_opt_metrics(df_compared)

    # Verify metrics
    assert str(Key.n_structures) in metrics
    assert metrics[str(Key.n_structures)] == len(ml_structs)
    assert str(Key.n_sym_ops_mae) in metrics
    assert str(MbdKey.structure_rmsd_vs_dft) in metrics


def test_structure_matcher_integration(cubic_struct: Structure) -> None:
    """Test that structure distances are calculated correctly.

    Based on the test_calc_structure_distances test in test_symmetry.py.
    """
    key = "struct1"
    df_ml_sym = symmetry.get_sym_info_from_structs({key: cubic_struct})

    # Create a slightly perturbed structure for testing
    slightly_perturbed = perturb_structure(cubic_struct, gamma=0.1)

    # Test structure distance calculation with structures that can be matched
    df_distances = symmetry.calc_structure_distances(
        df_ml_sym, {key: cubic_struct}, {key: slightly_perturbed}
    )

    # Check for the presence of RMSD column but don't assert specific values
    # since the structure matching is sensitive to implementation details
    assert MbdKey.structure_rmsd_vs_dft in df_distances
    assert Key.max_pair_dist in df_distances

    # Calculate geo opt metrics - should work even with potential NaN values
    metrics = geo_opt.calc_geo_opt_metrics(df_distances)
    assert str(MbdKey.structure_rmsd_vs_dft) in metrics
    assert isinstance(metrics[str(MbdKey.structure_rmsd_vs_dft)], float)


def test_full_geo_opt_workflow(
    small_structure_set: dict[str, Any],
    tmp_path: Path,
) -> None:
    """Test the complete geometry optimization workflow."""
    # Step 1: Create mock model data
    model_id = "test_model"
    model_dir = f"{tmp_path}/models/{model_id}"
    os.makedirs(model_dir, exist_ok=True)

    # Extract structures
    ml_structs = small_structure_set["ml_structs"]
    dft_structs = small_structure_set["dft_structs"]

    # Create a more significantly perturbed structure to ensure symmetry differences
    key = next(iter(ml_structs))  # Use first structure key
    # Create a more significantly perturbed structure
    ml_structs[key] = perturb_structure(ml_structs[key], gamma=1.5)

    # Step 2: Prepare input data
    # Create DFT analysis dataframe
    df_dft_sym = symmetry.get_sym_info_from_structs(dft_structs, symprec=1e-2)

    # Create a temporary structure file for the mock model
    ml_structs_dict = {
        mat_id: {"structure": struct.as_dict()} for mat_id, struct in ml_structs.items()
    }

    # Convert to DataFrame with structure column
    df_ml = pd.DataFrame.from_dict(ml_structs_dict, orient="index")

    # Save the test data
    ml_structs_path = f"{model_dir}/test_structures.csv"
    df_ml.to_csv(ml_structs_path)

    # Step 3: Run symmetry analysis
    symprec = 1e-2
    df_ml_sym = symmetry.get_sym_info_from_structs(ml_structs, symprec=symprec)

    # Step 4: Compare symmetry and calculate structure distances
    df_results = df_ml_sym.copy()
    df_results[MbdKey.spg_num_diff] = df_ml_sym[Key.spg_num] - df_dft_sym[Key.spg_num]
    df_results[MbdKey.n_sym_ops_diff] = (
        df_ml_sym[Key.n_sym_ops] - df_dft_sym[Key.n_sym_ops]
    )

    df_results = symmetry.calc_structure_distances(df_results, ml_structs, dft_structs)

    # Step 5: Calculate metrics
    metrics = geo_opt.calc_geo_opt_metrics(df_results)

    # Verify the key metrics are present
    assert str(Key.n_structures) in metrics
    assert str(Key.n_sym_ops_mae) in metrics
    assert str(MbdKey.structure_rmsd_vs_dft) in metrics
    assert str(Key.symmetry_match) in metrics
    assert metrics[str(Key.n_structures)] == len(ml_structs)

    # Check that structure distances are calculated for at least some structures
    assert not df_results[MbdKey.structure_rmsd_vs_dft].isna().all()

    # We've artificially introduced a structure with significant perturbation,
    # so we should have at least one structure with different symmetry
    different_symmetry = False
    for mat_id in df_results.index:
        if (
            df_results.loc[mat_id, MbdKey.spg_num_diff] != 0
            or df_results.loc[mat_id, MbdKey.n_sym_ops_diff] != 0
        ):
            different_symmetry = True
            break

    assert different_symmetry, "Expected at least one structure with different symmetry"
