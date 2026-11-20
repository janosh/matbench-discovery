"""Raw FC2/FC3 finite-displacement force sets of kappa runs, stored as HDF5.

Each material is one HDF5 group holding ``matscikit.AnharmonicPhonons.force_set()``
verbatim: the relaxed unit cell, the FC2/FC3 supercell and primitive matrices,
symprec and per order (``fc2``, plus ``fc3`` for materials whose FC3 stage ran) the
supercell geometry the atom indices refer to, the displacement distance, the
displacement datasets exactly as evaluated (FC3 pairs columnar, as phono3py's), the
forces on every displaced supercell (float32, eV/Å, fletcher32-checksummed) and how
much symmetrizing the fit changed it; the ``displacement_source`` attribute names the
code that chose the displacements. Datasets are replayed verbatim rather than
regenerated, so files stay loadable when a code chooses its displacements
differently. Forces are the source of truth: force constants follow from them by a
deterministic fit, but not vice versa. Float32 storage changes refitted force
constants by ~5e-8 (FC2) and ~2e-6 (FC3) of their largest element (EMT Cu3Au),
orders of magnitude below model errors.

Rebuild a material's fitted force constants with::

    phonons = load_phonons(path, material_id)  # matscikit.AnharmonicPhonons
    phonons.fc2, phonons.fc3_dense  # FC3 only for materials whose FC3 stage ran
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any
from urllib.parse import quote, unquote

import h5py
import matscikit
import numpy as np
from matscikit import AnharmonicPhonons

from matbench_discovery import atomic_output_path

# 3: matscikit.AnharmonicPhonons.force_set() layout, checked by from_force_set
# (scripts/migrate_force_sets_v2_to_v3.py converts schema-2 files)
FORCE_SET_SCHEMA_VERSION = 3
FORCE_DTYPE = np.float32


def _write_group(group: h5py.Group, mapping: Mapping[str, Any]) -> None:
    """Nested dicts become groups, None entries (no FC3, unsymmetrized fits) are left
    out, forces are float32 with shuffle + gzip and a fletcher32 checksum h5py verifies
    on every read.
    """
    for key, value in mapping.items():
        if isinstance(value, Mapping):
            _write_group(group.create_group(key), value)
        elif key in ("forces", "pair_forces"):
            group.create_dataset(
                key,
                data=np.asarray(value, dtype=FORCE_DTYPE),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
                fletcher32=True,
            )
        elif value is not None:
            group.create_dataset(key, data=value)


def _read_group(group: h5py.Group) -> dict[str, Any]:
    """Inverse of ``_write_group`` (without the left-out None entries)."""
    return {
        key: _read_group(item) if isinstance(item, h5py.Group) else item[()]
        for key, item in group.items()
    }


def write_force_set(
    group: h5py.Group,
    force_set: Mapping[str, Any],
    *,
    displacement_source: str = f"matscikit {matscikit.__version__}",
) -> None:
    """Write a ``force_set()`` dict into an empty HDF5 group, recording the code that
    chose its displacements.
    """
    group.attrs["displacement_source"] = displacement_source
    _write_group(group, force_set)


def read_force_set(group: h5py.Group) -> dict[str, Any]:
    """Read one material's force set as ``AnharmonicPhonons.from_force_set`` takes it
    (forces as float32 arrays, ``fc3`` None when absent) plus its
    ``displacement_source``.
    """
    return {
        "fc3": None,
        **_read_group(group),
        "displacement_source": str(group.attrs["displacement_source"]),
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
    """Return the material IDs stored in a merged force-set file, sorted (HDF5 lists
    groups alphabetically, not in insertion order).
    """
    with h5py.File(path, "r") as h5_file:
        _check_schema(h5_file, path)
        return sorted(unquote(name) for name in h5_file)


def load_force_set(path: str, material_id: str) -> dict[str, Any]:
    """Read one material's force set from a merged force-set file."""
    with h5py.File(path, "r") as h5_file:
        _check_schema(h5_file, path)
        group_name = material_group_name(material_id)
        if group_name not in h5_file:
            raise KeyError(f"{material_id!r} not in {path}")
        return read_force_set(h5_file[group_name])


def load_phonons(path: str, material_id: str) -> AnharmonicPhonons:
    """Rebuild a material's phonons and refit FC2 (and FC3 when stored) from forces.

    The stored displacement datasets are replayed verbatim, so the forces are fitted
    to the exact supercells they were computed on; matscikit rejects geometry,
    displacement-length and mass mismatches. FC3 data is absent for materials whose
    conductivity stage was gated off.
    """
    return AnharmonicPhonons.from_force_set(load_force_set(path, material_id))
