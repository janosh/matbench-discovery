"""Raw FC2/FC3 finite-displacement force sets of kappa runs, stored as HDF5.

Each material is one HDF5 group holding everything needed to rebuild its force
constants with phono3py: the relaxed unit cell, the FC2/FC3 supercell and primitive
matrices, the phono3py displacement datasets (JSON, forces stripped) and the forces
on every displaced supercell (float32, eV/Å). Forces are the source of truth:
force constants follow from them by a deterministic phono3py fit, but not vice
versa. Float32 storage changes refitted force constants by ~5e-8 (FC2) and ~1e-6
(FC3) of their largest element (EMT Al), orders of magnitude below model errors.

Rebuild a material's force constants with::

    ph3 = load_phono3py(path, material_id)
    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.produce_fc3(symmetrize_fc3r=True)  # materials whose FC3 stage ran
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING, Any
from urllib.parse import quote, unquote

import h5py
import numpy as np

from matbench_discovery import atomic_output_path

if TYPE_CHECKING:
    from phono3py import Phono3py

FORCE_SET_SCHEMA_VERSION = 1
FORCE_DTYPE = np.float32
_CELL_KEYS = ("cell", "scaled_positions", "numbers", "masses")
_MATRIX_KEYS = ("supercell_matrix", "phonon_supercell_matrix", "primitive_matrix")


def _strip_forces(dataset: object) -> object:
    """Return a copy of a phono3py displacement dataset without its forces.

    phono3py stores the assigned forces inside its datasets as well, which would
    duplicate the separately stored force arrays.
    """
    if isinstance(dataset, Mapping):
        return {
            str(key): _strip_forces(value)
            for key, value in dataset.items()
            if key != "forces"
        }
    if isinstance(dataset, list | tuple):
        return [_strip_forces(item) for item in dataset]
    if isinstance(dataset, np.ndarray):
        return dataset.tolist()
    if isinstance(dataset, np.generic):
        return dataset.item()
    return dataset


def force_set_from_phono3py(ph3: Phono3py, fc2_forces: np.ndarray) -> dict[str, Any]:
    """Collect one material's FC2 force set and the cell it was computed on.

    Args:
        ph3: Phono3py object whose FC2 displacements were evaluated.
        fc2_forces: (n_fc2_displacements, n_fc2_supercell_atoms, 3) forces (eV/Å).

    Returns:
        dict[str, Any]: force-set data; add_fc3_forces extends it with FC3 forces.
    """
    unitcell = ph3.unitcell
    return {
        "unitcell": {
            "cell": np.asarray(unitcell.cell),
            "scaled_positions": np.asarray(unitcell.scaled_positions),
            "numbers": np.asarray(unitcell.numbers),
            "masses": np.asarray(unitcell.masses),
        },
        "supercell_matrix": np.asarray(ph3.supercell_matrix),
        "phonon_supercell_matrix": np.asarray(ph3.phonon_supercell_matrix),
        "primitive_matrix": np.asarray(ph3.primitive_matrix),
        "symprec": float(ph3.symmetry.tolerance),
        "fc2_dataset": _strip_forces(ph3.phonon_dataset),
        "fc2_forces": np.asarray(fc2_forces),
        "fc3_dataset": None,
        "fc3_forces": None,
    }


def add_fc3_forces(
    force_set: Mapping[str, Any], ph3: Phono3py, fc3_forces: np.ndarray
) -> dict[str, Any]:
    """Return force_set extended with the FC3 displacement dataset and forces."""
    return {
        **force_set,
        "fc3_dataset": _strip_forces(ph3.dataset),
        "fc3_forces": np.asarray(fc3_forces),
    }


def _write_forces(group: h5py.Group, name: str, forces: np.ndarray) -> None:
    """Write a float32 force array with shuffle + gzip compression."""
    forces = np.asarray(forces, dtype=float)
    if forces.ndim != 3 or forces.shape[2] != 3 or not np.isfinite(forces).all():
        raise ValueError(f"{name} must be finite (n_disp, n_atoms, 3), {forces.shape=}")
    group.create_dataset(
        name,
        data=forces.astype(FORCE_DTYPE),
        compression="gzip",
        compression_opts=4,
        shuffle=True,
    )


def write_force_set(group: h5py.Group, force_set: Mapping[str, Any]) -> None:
    """Write one material's force set into an empty HDF5 group."""
    cell_group = group.create_group("unitcell")
    for key in _CELL_KEYS:
        cell_group.create_dataset(key, data=np.asarray(force_set["unitcell"][key]))
    for key in _MATRIX_KEYS:
        group.create_dataset(key, data=np.asarray(force_set[key], dtype=float))
    group.attrs["symprec"] = float(force_set["symprec"])
    for order in ("fc2", "fc3"):
        dataset, forces = force_set[f"{order}_dataset"], force_set[f"{order}_forces"]
        if (dataset is None) != (forces is None):
            raise ValueError(f"{order} dataset and forces must be given together")
        if forces is None:
            continue
        group.attrs[f"{order}_dataset"] = json.dumps(
            _strip_forces(dataset), sort_keys=True, allow_nan=False
        )
        _write_forces(group, f"{order}_forces", forces)
    if "fc2_forces" not in group:
        raise ValueError("A force set requires FC2 forces")


def read_force_set(group: h5py.Group) -> dict[str, Any]:
    """Read one material's force set (forces as float32 arrays)."""
    return {
        "unitcell": {key: group["unitcell"][key][()] for key in _CELL_KEYS},
        **{key: group[key][()] for key in _MATRIX_KEYS},
        "symprec": float(group.attrs["symprec"]),
        **{
            key: (
                json.loads(group.attrs[f"{order}_dataset"])
                if key.endswith("dataset")
                else group[f"{order}_forces"][()]
            )
            if f"{order}_forces" in group
            else None
            for order in ("fc2", "fc3")
            for key in (f"{order}_dataset", f"{order}_forces")
        },
    }


def _stamp_file(h5_file: h5py.File) -> None:
    """Record the schema version and units on a force-set file."""
    h5_file.attrs["schema_version"] = FORCE_SET_SCHEMA_VERSION
    h5_file.attrs["force_unit"] = "eV/Å"
    h5_file.attrs["length_unit"] = "Å"


def _check_schema(h5_file: h5py.File, path: str) -> None:
    """Reject force-set files from another schema version."""
    version = h5_file.attrs.get("schema_version")
    if version != FORCE_SET_SCHEMA_VERSION:
        raise ValueError(
            f"{path} has force-set schema {version}, "
            f"expected {FORCE_SET_SCHEMA_VERSION}"
        )


def material_group_name(material_id: str) -> str:
    """Return the HDF5 group name of a material (IDs may contain '/')."""
    return quote(material_id, safe="-_.")


def write_force_set_file(path: str, force_set: Mapping[str, Any]) -> None:
    """Atomically write a single-material force-set file (shard sidecar)."""
    with (
        atomic_output_path(path) as temporary_path,
        h5py.File(temporary_path, "w") as h5_file,
    ):
        _stamp_file(h5_file)
        write_force_set(h5_file.create_group("material"), force_set)


def merge_force_set_files(path: str, sources: Iterable[tuple[str, str]]) -> None:
    """Atomically merge (material_id, single-material file) pairs into one file."""
    with (
        atomic_output_path(path) as temporary_path,
        h5py.File(temporary_path, "w") as merged,
    ):
        _stamp_file(merged)
        for material_id, source_path in sources:
            with h5py.File(source_path, "r") as source:
                _check_schema(source, source_path)
                source.copy(
                    source["material"], merged, material_group_name(material_id)
                )


def material_ids(path: str) -> list[str]:
    """Return the material IDs stored in a merged force-set file."""
    with h5py.File(path, "r") as h5_file:
        _check_schema(h5_file, path)
        return [unquote(name) for name in h5_file]


def load_force_set(path: str, material_id: str) -> dict[str, Any]:
    """Read one material's force set from a merged force-set file."""
    with h5py.File(path, "r") as h5_file:
        _check_schema(h5_file, path)
        group_name = material_group_name(material_id)
        if group_name not in h5_file:
            raise KeyError(f"{material_id!r} not in {path}")
        return read_force_set(h5_file[group_name])


def load_phono3py(path: str, material_id: str) -> Phono3py:
    """Rebuild a material's Phono3py object with its stored datasets and forces.

    Call produce_fc2 / produce_fc3 on the result to refit the force constants.
    FC3 data is absent for materials whose conductivity stage was gated off.
    """
    from phono3py import Phono3py
    from phonopy.structure.atoms import PhonopyAtoms

    force_set = load_force_set(path, material_id)
    cell = force_set["unitcell"]
    ph3 = Phono3py(
        PhonopyAtoms(
            cell=cell["cell"],
            scaled_positions=cell["scaled_positions"],
            numbers=cell["numbers"],
            masses=cell["masses"],
        ),
        supercell_matrix=force_set["supercell_matrix"],
        phonon_supercell_matrix=force_set["phonon_supercell_matrix"],
        primitive_matrix=force_set["primitive_matrix"],
        symprec=force_set["symprec"],
    )
    ph3.phonon_dataset = force_set["fc2_dataset"]
    ph3.phonon_forces = np.asarray(force_set["fc2_forces"], dtype=float)
    if force_set["fc3_forces"] is not None:
        ph3.dataset = force_set["fc3_dataset"]
        ph3.forces = np.asarray(force_set["fc3_forces"], dtype=float)
    return ph3
