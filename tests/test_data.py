import os
import sys
import zipfile
from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from ferrox import Structure
from pymatviz.enums import Key
from ruamel.yaml.comments import CommentedMap

from matbench_discovery.data import (
    as_dict_handler,
    df_wbm,
    glob_to_df,
    load_df_wbm_with_preds,
    round_trip_yaml,
    structures_from_zip,
    structures_to_zip,
    update_yaml_file,
)
from matbench_discovery.enums import MbdKey, Model, TestSubset

struct_h2o = Structure(
    lattice=[5, 5, 5, 90, 90, 90],
    species=["H", "H", "O"],
    xyz=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],
)
struct_co2 = Structure(
    lattice=[6, 6, 6, 90, 90, 90],
    species=["C", "O", "O"],
    xyz=[[0, 0, 0], [0, 0, 1], [0, 0, -1]],
)

dummy_structures: dict[str, Structure] = {
    "structure1": struct_h2o,
    "structure2": struct_co2,
}


def test_as_dict_handler() -> None:
    """Test as_dict_handler serialization."""

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
    """Test WBM summary dataframe shape and columns."""
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
    """Test glob_to_df combines files matching pattern."""
    os.makedirs(tmp_path, exist_ok=True)
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


def test_structures_zip_round_trip(tmp_path: Path) -> None:
    """Test writing and reading structures to/from ZIP archive."""
    zip_path = tmp_path / "test_structures.zip"
    structures_to_zip(dummy_structures, zip_path)

    read_structs = structures_from_zip(zip_path)

    assert len(read_structs) == len(dummy_structures)

    for mat_id, original in dummy_structures.items():
        assert mat_id in read_structs, f"{mat_id=} not in read_structs"
        read = read_structs[mat_id]
        assert original.formula == read.formula
        assert np.allclose(original.cart_coords, read.cart_coords, atol=1e-6)
        assert len(original) == len(read)


def test_structures_from_zip_with_file_filter(tmp_path: Path) -> None:
    """Test file_filter parameter on structures_from_zip."""
    zip_path = tmp_path / "test_structures.zip"
    structures_to_zip(dummy_structures, zip_path)
    read_structs = structures_from_zip(
        zip_path, file_filter=lambda name, _idx: "1" in name
    )
    assert len(read_structs) == 1
    assert "structure1" in read_structs


def test_structures_from_zip_empty_file(tmp_path: Path) -> None:
    """Test reading from an empty ZIP file."""
    empty_zip = tmp_path / "empty.zip"
    with zipfile.ZipFile(empty_zip, mode="w"):
        pass
    read_structs = structures_from_zip(empty_zip)
    assert len(read_structs) == 0


def test_structures_from_zip_invalid_file(tmp_path: Path) -> None:
    """Test reading from an invalid ZIP file raises BadZipFile."""
    invalid_zip = tmp_path / "invalid.zip"
    with open(invalid_zip, mode="w") as file:
        file.write("This is not a zip file")
    with pytest.raises(zipfile.BadZipFile):
        structures_from_zip(invalid_zip)


def test_structures_from_zip_with_limit(tmp_path: Path) -> None:
    """Test limit parameter on structures_from_zip."""
    zip_path = tmp_path / "test_structures.zip"
    structures_to_zip(dummy_structures, zip_path)

    # Test with limit=1
    read_structs = structures_from_zip(zip_path, limit=1)
    assert len(read_structs) == 1

    # Test with limit=2 (should read all structures as there are only 2)
    read_structs = structures_from_zip(zip_path, limit=2)
    assert len(read_structs) == 2

    # Test with limit=None (default behavior, should read all structures)
    read_structs = structures_from_zip(zip_path, limit=None)
    assert len(read_structs) == 2

    # Test with limit greater than the number of structures
    read_structs = structures_from_zip(zip_path, limit=10)
    assert len(read_structs) == 2


@pytest.mark.parametrize(
    "slice_limit, expected_count",
    [
        (slice(1, 3), 2),
        (slice(None, 2), 2),
        (slice(3, None), 2),
        (slice(0, None, 2), 3),
        (slice(None, None, -1), 5),
        (slice(5, 5), 0),
        (slice(10, 20), 0),
    ],
)
def test_structures_from_zip_with_slice_limit(
    tmp_path: Path,
    slice_limit: slice,
    expected_count: int,
) -> None:
    """Test structures_from_zip with slice objects for the limit parameter."""
    more_structures: dict[str, Structure] = {
        "slice_test_0": Structure(
            lattice=[4, 4, 4, 90, 90, 90],
            species=["H", "H"],
            xyz=[[0, 0, 0], [0, 0, 1]],
        ),
        "slice_test_1": Structure(
            lattice=[4, 4, 4, 90, 90, 90],
            species=["N", "N"],
            xyz=[[0, 0, 0], [0, 0, 1.1]],
        ),
        "slice_test_2": Structure(
            lattice=[4, 4, 4, 90, 90, 90],
            species=["O", "O"],
            xyz=[[0, 0, 0], [0, 0, 1.2]],
        ),
        "slice_test_3": Structure(
            lattice=[4, 4, 4, 90, 90, 90],
            species=["F", "F"],
            xyz=[[0, 0, 0], [0, 0, 1.4]],
        ),
        "slice_test_4": Structure(
            lattice=[5, 5, 5, 90, 90, 90],
            species=["Cl", "Cl"],
            xyz=[[0, 0, 0], [0, 0, 2.0]],
        ),
    }

    zip_path = tmp_path / "slice_test_structures.zip"
    structures_to_zip(more_structures, zip_path)

    read_structs = structures_from_zip(zip_path, limit=slice_limit)
    assert len(read_structs) == expected_count


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
    """Test loading WBM predictions with model filtering."""
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
    """Test max_error_threshold filtering in load_df_wbm_with_preds."""
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
        # Make glob return a non-empty list to skip the mock data loading path
        patch("matbench_discovery.data.glob", return_value=["dummy_file.csv"]),
        # Patch only the specific read_csv call in glob_to_df that loads predictions,
        # not the one that loads mock data
        patch("pandas.read_csv", return_value=df_float),
        pytest.raises(ValueError, match=r"pred_col.*not found in"),
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


def test_update_yaml_file(tmp_path: Path) -> None:
    """Test updating YAML files at specific paths."""
    test_file = f"{tmp_path}/test.yml"

    # Test case 1: Basic update at root level
    initial_data = {"metrics": {"discovery": {"mae": 0.1}}}
    with open(test_file, mode="w") as file:
        round_trip_yaml.dump(initial_data, file)

    updated_yaml = update_yaml_file(
        test_file, "metrics.discovery", {"mae": 0.2, "rmse": 0.3}
    )
    assert updated_yaml["metrics"]["discovery"] == {"mae": 0.2, "rmse": 0.3}

    # Test case 2: Create new nested path
    updated_yaml = update_yaml_file(test_file, "metrics.new.nested.path", {"value": 42})
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

    updated_yaml = update_yaml_file(
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
    updated_yaml = update_yaml_file(test_file, "new.path", commented_data)

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
        update_yaml_file("non-existent.yml", "path", {"data": 1})

    # Test bad paths
    for path in ("metrics..discovery", "metrics..", "metrics.discovery..", "."):
        with pytest.raises(ValueError, match="Invalid dotted_path="):
            update_yaml_file(test_file, path, {"data": 1})
