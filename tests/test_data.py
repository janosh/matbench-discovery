import os
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from pymatgen.core import Lattice, Structure
from pymatviz.enums import Key

from matbench_discovery.data import (
    as_dict_handler,
    ase_atoms_from_zip,
    ase_atoms_to_zip,
    df_wbm,
    glob_to_df,
)

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
def test_glob_to_df(pattern: str, tmp_path: Path, df_mixed: pd.DataFrame) -> None:
    os.makedirs(f"{tmp_path}", exist_ok=True)
    df_mixed.to_csv(f"{tmp_path}/dummy_df.csv", index=False)
    df_mixed.to_json(f"{tmp_path}/dummy_df.json")

    df_out = glob_to_df(f"{tmp_path}/{pattern}")
    assert df_out.shape == df_mixed.shape
    assert list(df_out) == list(df_mixed)

    with pytest.raises(FileNotFoundError):
        glob_to_df("foo")


def test_atoms_zip_round_trip(tmp_path: Path) -> None:
    # Write atoms to a temporary ZIP file
    zip_path = tmp_path / "test_structures.zip"
    ase_atoms_to_zip(dummy_atoms, zip_path)

    # Read atoms back from the ZIP file
    read_atoms = ase_atoms_from_zip(zip_path)

    # Check that we got the same number of Atoms objects back
    assert len(read_atoms) == len(dummy_atoms)

    for original, read in zip(dummy_atoms, read_atoms):
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
    with zipfile.ZipFile(empty_zip, "w"):
        pass
    read_atoms = ase_atoms_from_zip(empty_zip)
    assert len(read_atoms) == 0


def test_ase_atoms_from_zip_invalid_file(tmp_path: Path) -> None:
    invalid_zip = tmp_path / "invalid.zip"
    with open(invalid_zip, "w") as file:
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
