"""Convert schema-2 kappa force-set files to schema 3 without any model call.

Schema 2 stored matscikit-order (atom, vector) datasets as JSON attributes; schema 3 is
``matscikit.AnharmonicPhonons.force_set()`` written group by group (FC3 pairs columnar).
The float32 forces are copied bit for bit, so a migrated file refits exactly the force
constants the schema-2 loader fitted: symmetrized, recording each order's measured
``symmetrization`` (which ``from_force_set`` replays). Every material is reloaded with
``force_sets.load_phonons`` before the output is kept, so matscikit's geometry,
displacement-length and mass checks pass on all of them.

    uv run --with-editable . scripts/migrate_force_sets_v2_to_v3.py \
        2026-11-09-phonons-kappa-103-forces.h5 migrated-forces.h5
"""

import argparse
import json
from typing import Any

import h5py
import numpy as np
from matscikit import AnharmonicPhonons

from matbench_discovery import atomic_output_path
from matbench_discovery.phonons import force_sets


def _singles(dataset: dict[str, Any], displacements: list[Any]) -> dict[str, Any]:
    """Schema-3 single displacements from a schema-2 dataset's supercell geometry."""
    return {
        "lattice": np.asarray(dataset["supercell"]["lattice"], dtype=float),
        "frac_coords": np.asarray(dataset["supercell"]["frac_coords"], dtype=float),
        "atoms": np.array([atom for atom, _vector in displacements], dtype=np.int64),
        "vectors": np.array([vector for _atom, vector in displacements], dtype=float),
    }


def migrate_material(group: h5py.Group) -> dict[str, Any]:
    """The schema-3 force set of one schema-2 material group."""
    cell = group["unitcell"]
    distance = float(group.attrs["displacement_distance"])
    fc2_dataset = json.loads(group.attrs["fc2_dataset"])
    force_set: dict[str, Any] = {
        "unitcell": {
            "lattice": cell["cell"][()],
            "frac_coords": cell["scaled_positions"][()],
            "numbers": cell["numbers"][()],
            "masses": cell["masses"][()],
        },
        # schema 2 stored M with supercell lattice = M @ unit-cell lattice, the
        # transpose of the constructor's (phonopy's) supercell matrix
        **{
            key: np.rint(group[key][()]).astype(np.int64).T
            for key in ("supercell_matrix", "phonon_supercell_matrix")
        },
        "primitive_matrix": group["primitive_matrix"][()],
        "symprec": float(group.attrs["symprec"]),
        "fc2": {
            **_singles(fc2_dataset, fc2_dataset["displacements"]),
            "distance": distance,
            "forces": group["fc2_forces"][()],
        },
        "fc3": None,
    }
    if "fc3_forces" in group:
        fc3_dataset = json.loads(group.attrs["fc3_dataset"])
        configs = fc3_dataset["displacements"]
        singles = [config[0] for config in configs if len(config) == 1]
        pairs = [config for config in configs if len(config) == 2]
        forces = group["fc3_forces"][()]
        # forces are split positionally: all singles first, then all pairs
        lengths = [len(config) for config in configs]
        expected = [1] * len(singles) + [2] * len(pairs)
        if lengths != expected or len(forces) != len(configs):
            raise ValueError(
                f"{group.name}: expected {len(singles)} single then {len(pairs)} pair "
                f"configs matching {len(forces)} force rows, got lengths {lengths}"
            )
        force_set["fc3"] = {
            **_singles(fc3_dataset, singles),
            "distance": distance,
            "forces": forces[: len(singles)],
            "pair_parents": np.array(
                [singles.index(first) for first, _second in pairs], dtype=np.int64
            ),
            "pair_atoms": np.array(
                [second[0] for _first, second in pairs], dtype=np.int64
            ),
            "pair_vectors": np.array(
                [second[1] for _first, second in pairs], dtype=float
            ).reshape(-1, 3),
            "pair_forces": forces[len(singles) :],
        }
    return _record_symmetrization(force_set)


def _record_symmetrization(force_set: dict[str, Any]) -> dict[str, Any]:
    """Record the symmetrization of each stored order: schema 2 always fitted
    symmetrized, so refit with a finite marker (any finite value means symmetrized)
    and store the measured change in its place.
    """
    orders = [order for order in ("fc2", "fc3") if force_set[order] is not None]
    for order in orders:
        force_set[order]["symmetrization"] = 0.0
    measured = AnharmonicPhonons.from_force_set(force_set).symmetrization
    for order, value in zip(("fc2", "fc3"), measured, strict=True):
        if order in orders:
            force_set[order]["symmetrization"] = value
    return force_set


def migrate(source: str, target: str) -> None:
    """Atomically write the schema-3 conversion of schema-2 file ``source``, checking
    that every material loads.
    """
    with (
        h5py.File(source, "r") as old,
        atomic_output_path(target) as temporary_path,
    ):
        if (version := old.attrs.get("schema_version")) != 2:
            raise ValueError(f"{source} has force-set schema {version}, expected 2")
        with h5py.File(temporary_path, "w") as new:
            force_sets._stamp_file(new)  # noqa: SLF001
            for name, group in old.items():
                force_sets.write_force_set(
                    new.create_group(name),
                    migrate_material(group),
                    displacement_source=str(group.attrs["displacement_source"]),
                )
        for material_id in force_sets.material_ids(temporary_path):
            force_sets.load_phonons(temporary_path, material_id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", help="schema-2 force-set HDF5 file")
    parser.add_argument("target", help="schema-3 output path")
    args = parser.parse_args()
    migrate(args.source, args.target)
