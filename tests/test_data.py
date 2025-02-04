import os
import sys
import zipfile
from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import requests
from ase import Atoms
from pymatgen.core import Lattice, Structure
from pymatviz.enums import Key

from matbench_discovery import DATA_DIR
from matbench_discovery.data import (
    DataFiles,
    Files,
    Model,
    as_dict_handler,
    ase_atoms_from_zip,
    ase_atoms_to_zip,
    df_wbm,
    glob_to_df,
    load_df_wbm_with_preds,
)
from matbench_discovery.enums import MbdKey, TestSubset

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
    read_atoms = ase_atoms_from_zip(zip_path, file_filter=lambda name: "1" in name)
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


def test_files() -> None:
    """Test error handling in Files enum."""

    assert Files.base_dir == DATA_DIR

    # Test custom base_dir
    class SubFiles(Files, base_dir="foo"):
        pass

    assert SubFiles.base_dir == "foo"

    # Test invalid label lookup
    label = "invalid-label"
    with pytest.raises(ValueError, match=f"{label=} not found in Files"):
        Files.from_label(label)


def test_data_files() -> None:
    """Test DataFiles enum functionality."""
    # Test that paths are constructed correctly
    assert str(DataFiles.mp_energies) == f"{DATA_DIR}/mp/2023-01-10-mp-energies.csv.gz"
    assert repr(DataFiles.mp_energies) == "DataFiles.mp_energies"
    assert DataFiles.mp_energies.name == "mp_energies"
    assert (
        DataFiles.mp_energies.url == "https://figshare.com/ndownloader/files/49083124"
    )
    assert DataFiles.mp_energies.rel_path == "mp/2023-01-10-mp-energies.csv.gz"

    # Test that multiple files exist and have correct attributes
    assert DataFiles.wbm_summary.rel_path == "wbm/2023-12-13-wbm-summary.csv.gz"
    assert DataFiles.wbm_summary.path == f"{DATA_DIR}/wbm/2023-12-13-wbm-summary.csv.gz"
    assert (
        DataFiles.wbm_summary.url == "https://figshare.com/ndownloader/files/44225498"
    )


def test_model() -> None:
    """Test Model enum functionality."""
    # Test basic model attributes
    assert Model.alignn.name == "alignn"
    assert Model.alignn.rel_path == "alignn/alignn.yml"
    assert Model.alignn.url is None
    assert Model.alignn.label == "ALIGNN"

    # Test metadata property
    metadata = Model.alignn.metadata
    assert isinstance(metadata, dict)

    # Test yaml_path property
    assert Model.alignn.yaml_path.endswith("alignn/alignn.yml")
    assert Model.grace2l_r6.kappa_103_path.endswith(
        "kappa-FIRE-fmax=1e-4-symprec=1e-5.json.gz"
    )

    # Test error handling for missing paths
    with pytest.raises(
        AttributeError,
        match="'Model' object has no attribute 'bad_path'",
    ):
        _ = Model.alignn.bad_path

    # Test Model metrics property
    metrics = Model.alignn.metrics
    assert isinstance(metrics, dict)


@pytest.mark.parametrize("data_file", DataFiles)
def test_data_files_urls(data_file: DataFiles) -> None:
    """Test that each URL in data-files.yml is a valid Figshare download URL."""

    name, url = data_file.name, data_file.url
    # check that URL is a figshare download
    assert "figshare.com/ndownloader/files/" in url, (
        f"URL for {name} is not a Figshare download URL: {url}"
    )

    # check that the URL is valid by sending a head request
    response = requests.head(url, allow_redirects=True, timeout=5)
    assert response.status_code in {200, 403}, f"Invalid URL for {name}: {url}"


def test_download_file(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    """Test download_file function."""

    from matbench_discovery.data import download_file

    url = "https://example.com/test.txt"
    test_content = b"test content"
    dest_path = tmp_path / "test.txt"

    # Mock successful request
    mock_response = requests.Response()
    mock_response.status_code = 200
    mock_response._content = test_content  # noqa: SLF001

    with patch("requests.get", return_value=mock_response):
        download_file(str(dest_path), url)
        assert dest_path.read_bytes() == test_content

    # Mock failed request
    mock_response = requests.Response()
    mock_response.status_code = 404
    mock_response._content = b"Not found"  # noqa: SLF001

    with patch("requests.get", return_value=mock_response):
        download_file(str(dest_path), url)  # Should print error but not raise

    stdout, stderr = capsys.readouterr()
    assert f"Error downloading {url=}" in stdout
    assert stderr == ""


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
