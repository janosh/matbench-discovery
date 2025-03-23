import os
import sys
import zipfile
from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from pymatgen.core import Lattice, Structure
from pymatviz.enums import Key
from ruamel.yaml.comments import CommentedMap

from matbench_discovery.data import (
    as_dict_handler,
    ase_atoms_from_zip,
    ase_atoms_to_zip,
    df_wbm,
    glob_to_df,
    load_df_wbm_with_preds,
    round_trip_yaml,
    update_yaml_at_path,
)
from matbench_discovery.enums import MbdKey, Model, TestSubset

structure = Structure(
    lattice=Lattice.cubic(5),
    species=("Fe", "O"),
    coords=((0, 0, 0), (0.5, 0.5, 0.5)),
)

atoms1 = Atoms("H2O", positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]], cell=[5, 5, 5])
atoms1.info[Key.mat_id] = "structure1"
atoms1.info["formation_energy"] = 1.23
atoms1.info["forces"] = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
atoms1.info["dict_data"] = {"a": 1, "b": 2}
atoms1.info["list_data"] = [1, 2, 3]

atoms2 = Atoms("CO2", positions=[[0, 0, 0], [0, 0, 1], [0, 0, -1]], cell=[6, 6, 6])
atoms2.info[Key.mat_id] = "structure2"
atoms2.info["bandgap"] = 2.34
atoms2.info["magmoms"] = np.array([0.1, -0.1, 0.1])

dummy_atoms = [atoms1, atoms2]


def test_as_dict_handler() -> None:
    class C:
        @staticmethod
        def as_dict() -> dict[str, Any]:
            return {"foo": "bar"}

    assert as_dict_handler(C()) == {"foo": "bar"}
    assert as_dict_handler(1) is None
    assert as_dict_handler("foo") is None
    assert as_dict_handler([1, 2, 3]) is None
    assert as_dict_handler({"foo": "bar"}) is None


def test_df_wbm() -> None:
    assert df_wbm.shape == (256_963, 18)
    assert df_wbm.index.name == Key.mat_id
    assert set(df_wbm) > {Key.formula, Key.mat_id, Key.bandgap_pbe}

    for col in (MbdKey.e_form_dft, MbdKey.each_true):
        assert col in df_wbm, f"{col=} not in {list(df_wbm)=}"


@pytest.mark.parametrize("pattern", ["*df.csv", "*df.json"])
def test_glob_to_df(
    pattern: str,
    tmp_path: Path,
    df_mixed: pd.DataFrame,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    os.makedirs(f"{tmp_path}", exist_ok=True)
    df_mixed.to_csv(f"{tmp_path}/dummy_df.csv", index=False)
    df_mixed.to_json(f"{tmp_path}/dummy_df.json")

    df_out = glob_to_df(f"{tmp_path}/{pattern}")
    assert df_out.shape == df_mixed.shape
    assert list(df_out) == list(df_mixed)

    with pytest.raises(ValueError, match="Unsupported file extension in pattern='foo'"):
        glob_to_df("foo")

    # Mock sys.modules without pytest to test file not found error
    mock_modules = dict(sys.modules)
    mock_modules.pop("pytest", None)  # remove pytest since glob_to_df returns mock data
    # if if finds pytest imported
    # also remove CI from os.environ
    monkeypatch.delenv("CI", raising=False)
    with (
        patch("sys.modules", mock_modules),
        pytest.raises(FileNotFoundError, match="No files matching glob pattern="),
    ):
        glob_to_df("foo.csv")


@pytest.mark.parametrize(
    "dummy_atoms",
    [dummy_atoms, {"atoms1": atoms1, "atoms2": atoms2}],
)
def test_atoms_zip_round_trip(
    tmp_path: Path, dummy_atoms: dict[str, Atoms] | list[Atoms]
) -> None:
    # Write atoms to a temporary ZIP file
    zip_path = tmp_path / "test_structures.zip"
    ase_atoms_to_zip(dummy_atoms, zip_path)

    # Read atoms back from the ZIP file
    read_atoms = ase_atoms_from_zip(zip_path)

    # Check that we got the same number of Atoms objects back
    assert len(read_atoms) == len(dummy_atoms)

    orig_atoms = dummy_atoms.values() if isinstance(dummy_atoms, dict) else dummy_atoms
    for original, read in zip(orig_atoms, read_atoms):
        # Check basic Atoms properties
        assert original.get_chemical_formula() == read.get_chemical_formula()
        assert np.allclose(original.get_positions(), read.get_positions())
        assert np.allclose(original.get_cell(), read.get_cell())
        assert np.all(original.pbc == read.pbc)

        # Check info dictionary
        for key, value in original.info.items():
            assert key in read.info, f"{key=} not in {list(read.info)}"
            read_val = read.info[key]
            if np.ndarray in {type(value), type(read_val)}:
                assert np.allclose(value, read_val), f"Mismatch in {key}"
            else:
                assert value == read_val, f"Mismatch in {key}"

        # Check that no extra keys were added
        assert set(original.info) == set(read.info)


def test_ase_atoms_from_zip_with_file_filter(tmp_path: Path) -> None:
    zip_path = tmp_path / "test_structures.zip"
    ase_atoms_to_zip(dummy_atoms, zip_path)
    read_atoms = ase_atoms_from_zip(
        zip_path, file_filter=lambda name, _idx: "1" in name
    )
    assert len(read_atoms) == 1
    assert read_atoms[0].get_chemical_formula() == "H2O"


def test_ase_atoms_from_zip_with_filename_to_info(tmp_path: Path) -> None:
    zip_path = tmp_path / "test_structures.zip"
    ase_atoms_to_zip(dummy_atoms, zip_path)
    read_atoms = ase_atoms_from_zip(zip_path, filename_to_info=True)
    assert all("filename" in atoms.info for atoms in read_atoms)
    assert read_atoms[0].info["filename"] == "structure1.extxyz"
    assert read_atoms[1].info["filename"] == "structure2.extxyz"


def test_ase_atoms_from_zip_empty_file(tmp_path: Path) -> None:
    empty_zip = tmp_path / "empty.zip"
    with zipfile.ZipFile(empty_zip, mode="w"):
        pass
    read_atoms = ase_atoms_from_zip(empty_zip)
    assert len(read_atoms) == 0


def test_ase_atoms_from_zip_invalid_file(tmp_path: Path) -> None:
    invalid_zip = tmp_path / "invalid.zip"
    with open(invalid_zip, mode="w") as file:
        file.write("This is not a zip file")
    with pytest.raises(zipfile.BadZipFile):
        ase_atoms_from_zip(invalid_zip)


def test_ase_atoms_from_zip_with_limit(tmp_path: Path) -> None:
    zip_path = tmp_path / "test_structures.zip"
    ase_atoms_to_zip(dummy_atoms, zip_path)

    # Test with limit=1
    read_atoms = ase_atoms_from_zip(zip_path, limit=1)
    assert len(read_atoms) == 1
    assert read_atoms[0].get_chemical_formula() == "H2O"

    # Test with limit=2 (should read all structures as there are only 2)
    read_atoms = ase_atoms_from_zip(zip_path, limit=2)
    assert len(read_atoms) == 2
    assert read_atoms[0].get_chemical_formula() == "H2O"
    assert read_atoms[1].get_chemical_formula() == "CO2"

    # Test with limit=None (default behavior, should read all structures)
    read_atoms = ase_atoms_from_zip(zip_path, limit=None)
    assert len(read_atoms) == 2

    # Test with limit greater than the number of structures
    read_atoms = ase_atoms_from_zip(zip_path, limit=10)
    assert len(read_atoms) == 2


@pytest.mark.parametrize(
    "slice_limit, expected_formulas, expected_count",
    [
        (slice(1, 3), ["N2", "O2"], 2),
        (slice(None, 2), ["H2", "N2"], 2),
        (slice(3, None), ["F2", "Cl2"], 2),
        (slice(0, None, 2), ["H2", "O2", "Cl2"], 3),
        (slice(None, None, -1), ["Cl2", "F2", "O2", "N2", "H2"], 5),
        (slice(5, 5), [], 0),
        (slice(10, 20), [], 0),
    ],
)
def test_ase_atoms_from_zip_with_slice_limit(
    tmp_path: Path,
    slice_limit: slice,
    expected_formulas: list[str],
    expected_count: int,
) -> None:
    """Test ase_atoms_from_zip with slice objects for the limit parameter."""
    # Create a zip file with more structures for better slice testing
    more_atoms = [
        Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]], cell=[4, 4, 4]),
        Atoms("N2", positions=[[0, 0, 0], [0, 0, 1.1]], cell=[4, 4, 4]),
        Atoms("O2", positions=[[0, 0, 0], [0, 0, 1.2]], cell=[4, 4, 4]),
        Atoms("F2", positions=[[0, 0, 0], [0, 0, 1.4]], cell=[4, 4, 4]),
        Atoms("Cl2", positions=[[0, 0, 0], [0, 0, 2.0]], cell=[5, 5, 5]),
    ]

    # Add material IDs to each structure
    for idx, atoms in enumerate(more_atoms):
        atoms.info[Key.mat_id] = f"slice_test_{idx}"

    zip_path = tmp_path / "slice_test_structures.zip"
    ase_atoms_to_zip(more_atoms, zip_path)

    # Test with the parametrized slice
    read_atoms = ase_atoms_from_zip(zip_path, limit=slice_limit)

    # Check the number of structures returned
    assert len(read_atoms) == expected_count

    # Check the chemical formulas of the returned structures
    for idx, formula in enumerate(expected_formulas):
        if idx < len(read_atoms):
            assert read_atoms[idx].get_chemical_formula() == formula


@pytest.mark.skipif(
    "CI" in os.environ,
    reason="CI uses mock data so don't check length against on-the-fly "
    "downloaded actual df_wbm",
)
@pytest.mark.parametrize("models", [[], ["wrenformer"]])
@pytest.mark.parametrize("max_error_threshold", [None, 5.0, 1.0])
def test_load_df_wbm_with_preds(
    models: list[str], max_error_threshold: float | None
) -> None:
    df_wbm_with_preds = load_df_wbm_with_preds(
        models=models, max_error_threshold=max_error_threshold
    )
    assert len(df_wbm_with_preds) == len(df_wbm)

    assert list(df_wbm_with_preds) == list(df_wbm) + [
        Model[model].label for model in models
    ]
    assert df_wbm_with_preds.index.name == Key.mat_id

    for model_name in models:
        model = Model[model_name]
        assert model.label in df_wbm_with_preds
        if max_error_threshold is not None:
            # Check if predictions exceeding the threshold are filtered out
            error = abs(
                df_wbm_with_preds[model.label] - df_wbm_with_preds[MbdKey.e_form_dft]
            )
            assert np.all(error[~error.isna()] <= max_error_threshold)
        else:
            # If no threshold is set, all predictions should be present
            assert df_wbm_with_preds[model.label].isna().sum() == 0


@pytest.mark.skipif(
    "CI" in os.environ, reason="CI uses mock data where other error thresholds apply"
)
def test_load_df_wbm_max_error_threshold() -> None:
    models: dict[str, int] = {  # map model to number of max allowed missing preds
        Model.mace_mp_0.label: 38  # before error is raised
    }
    df_no_thresh = load_df_wbm_with_preds(models=list(models))
    df_high_thresh = load_df_wbm_with_preds(models=list(models), max_error_threshold=10)
    df_low_thresh = load_df_wbm_with_preds(models=list(models), max_error_threshold=0.1)

    for model, n_missing in models.items():
        assert df_no_thresh[model].isna().sum() == n_missing
        assert df_high_thresh[model].isna().sum() <= df_no_thresh[model].isna().sum()
        assert df_high_thresh[model].isna().sum() <= df_low_thresh[model].isna().sum()


def test_load_df_wbm_with_preds_errors(df_float: pd.DataFrame) -> None:
    """Test error handling in load_df_wbm_with_preds function."""

    # Test invalid model name
    with pytest.raises(ValueError, match="expected subset of"):
        load_df_wbm_with_preds(models=["InvalidModel"])

    # Test negative error threshold
    with pytest.raises(
        ValueError, match="max_error_threshold must be a positive number"
    ):
        load_df_wbm_with_preds(max_error_threshold=-1)

    # Test pred_col not in predictions file
    with (
        patch("pandas.read_csv", return_value=df_float),
        pytest.raises(ValueError, match="pred_col.*not found in"),
    ):
        load_df_wbm_with_preds(models=["alignn"])


@pytest.mark.parametrize(
    "subset",
    ["unique_prototypes", TestSubset.uniq_protos, ["wbm-1-1", "wbm-1-2"], None],
)
def test_load_df_wbm_with_preds_subset(subset: Any) -> None:
    """Test subset handling in load_df_wbm_with_preds."""
    df_wbm = load_df_wbm_with_preds(subset=subset)
    assert isinstance(df_wbm, pd.DataFrame)


def test_update_yaml_at_path(tmp_path: Path) -> None:
    """Test updating YAML files at specific paths."""
    test_file = f"{tmp_path}/test.yml"

    # Test case 1: Basic update at root level
    initial_data = {"metrics": {"discovery": {"mae": 0.1}}}
    with open(test_file, mode="w") as file:
        round_trip_yaml.dump(initial_data, file)

    updated_yaml = update_yaml_at_path(
        test_file, "metrics.discovery", {"mae": 0.2, "rmse": 0.3}
    )
    assert updated_yaml["metrics"]["discovery"] == {"mae": 0.2, "rmse": 0.3}

    # Test case 2: Create new nested path
    updated_yaml = update_yaml_at_path(
        test_file, "metrics.new.nested.path", {"value": 42}
    )
    assert updated_yaml["metrics"]["new"]["nested"]["path"] == {"value": 42}

    # Test case 3: Update with comments
    yaml_with_comments = """
metrics:
  discovery:  # Discovery metrics
    mae: 0.1  # Mean absolute error
    rmse: 0.2  # Root mean squared error
"""
    with open(test_file, mode="w") as file:
        file.write(yaml_with_comments)

    updated_yaml = update_yaml_at_path(
        test_file, "metrics.discovery", {"mae": 0.3, "rmse": 0.4}
    )
    assert updated_yaml["metrics"]["discovery"] == {"mae": 0.3, "rmse": 0.4}

    # Verify comments are preserved in the file
    with open(test_file) as file:
        content = file.read()
    assert "discovery:  # Discovery metrics\n    mae: 0.3" in content, f"{content=}"

    # Test case 4: Update with CommentedMap
    commented_data = CommentedMap({"value": 1})
    commented_data.yaml_add_eol_comment("A comment", "value")
    updated_yaml = update_yaml_at_path(test_file, "new.path", commented_data)

    # Verify the data structure
    assert updated_yaml["new"]["path"]["value"] == 1
    # Verify comments in the file
    with open(test_file) as file:
        content = file.read()
    # check that old content is still there
    assert "discovery:  # Discovery metrics\n    mae: 0.3" in content, f"{content=}"
    # check new content was added
    assert "value: 1  # A comment" in content, f"{content=}"

    # Test case 5: Error cases
    with pytest.raises(FileNotFoundError):
        update_yaml_at_path("non-existent.yml", "path", {"data": 1})

    # Test bad paths
    for path in ("metrics..discovery", "metrics..", "metrics.discovery..", "."):
        with pytest.raises(ValueError, match="Invalid dotted_path="):
            update_yaml_at_path(test_file, path, {"data": 1})
