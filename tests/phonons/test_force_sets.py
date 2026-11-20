"""Tests for HDF5 FC2/FC3 force-set sidecars."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pytest
from matscikit import EMT, AnharmonicPhonons, Structure

from matbench_discovery.phonons import force_sets
from matbench_discovery.phonons import thermal_conductivity as ltc
from matbench_discovery.phonons.adapters import StandardKappaAdapter
from matbench_discovery.phonons.pipeline import KappaSettings


@pytest.fixture(scope="module")
def cu_calc() -> ltc.KappaCalculation:
    """EMT FC2/FC3 fits of primitive Cu, keeping their displacement forces."""
    lattice = 3.61 / 2 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    cu_atoms = Structure(lattice, ["Cu"], [[0, 0, 0]], fractional=True)
    # 3x3x3 FC3: in 2x2x2 every atom is its own inversion image, forcing FC3 = 0
    calc = ltc.init_phonons(
        cu_atoms,
        fc2_supercell=np.diag([2, 2, 2]),
        fc3_supercell=np.diag([3, 3, 3]),
        q_point_mesh=(4, 4, 4),
        displacement_distance=0.03,
    )
    StandardKappaAdapter().calculate_fc2(calc, EMT(), KappaSettings())
    StandardKappaAdapter().calculate_fc3(calc, EMT(), KappaSettings())
    return calc


def write_merged(tmp_path: Path, force_set: dict[str, Any], material_id: str) -> str:
    """Write force_set as a shard sidecar and merge it into a model file."""
    shard_path = f"{tmp_path}/shard.h5"
    force_sets.write_force_set_file(shard_path, force_set)
    merged_path = f"{tmp_path}/2026-11-09-phonons-kappa-103-forces.h5"
    force_sets.merge_force_set_files(merged_path, [(material_id, shard_path)])
    return merged_path


@pytest.mark.parametrize("with_fc3", [True, False])
def test_force_set_roundtrip_refits_force_constants(
    tmp_path: Path, cu_calc: ltc.KappaCalculation, with_fc3: bool
) -> None:
    """Merged HDF5 force sets rebuild FC2 (and FC3) to float32 force precision;
    FC3-gated materials store FC2 only.
    """
    reference = cu_calc.phonons
    force_set = reference.force_set()
    if not with_fc3:
        force_set["fc3"] = None
    merged_path = write_merged(tmp_path, force_set, "mp-30/x")
    assert force_sets.material_ids(merged_path) == ["mp-30/x"]
    stored = force_sets.load_force_set(merged_path, "mp-30/x")
    assert stored["fc2"]["forces"].dtype == np.float32
    assert stored["fc2"]["distance"] == 0.03
    assert stored["displacement_source"].startswith("matscikit ")
    assert (stored["fc3"] is None) is not with_fc3

    phonons = force_sets.load_phonons(merged_path, "mp-30/x")
    # float32 keeps ~7 significant digits of each force; measured refit errors are
    # 1.2e-8 (FC2) and 7.4e-7 (FC3) of each tensor's largest element, so 1e-5 leaves
    # >10x margin while still catching a lost or misassigned displacement
    pairs = [(phonons.fc2, reference.fc2)]
    if with_fc3:
        pairs.append((phonons.fc3_dense, reference.fc3_dense))
    else:
        assert phonons.fc3_forces is None
    for refit, fitted in pairs:
        expected = np.asarray(fitted)
        assert np.abs(expected).max() > 0.1  # FC3 must not vanish by symmetry
        np.testing.assert_allclose(
            refit, expected, rtol=0, atol=1e-5 * np.abs(expected).max()
        )


@pytest.mark.parametrize(
    ("mutate", "match"),
    [
        # a different cell makes stored atom indices refer to other atoms
        (
            lambda fs: fs["fc2"].update(lattice=1.01 * fs["fc2"]["lattice"]),
            "fc2 supercell lattice differs",
        ),
        # vectors no longer have the recorded length (unit or basis mix-up)
        (lambda fs: fs["fc2"].update(distance=0.01), "not the distance 0.01"),
        # FC3 rows refer to the FC3 supercell's atom order, checked separately
        (
            lambda fs: fs["fc3"].update(
                frac_coords=np.roll(fs["fc3"]["frac_coords"], 1, axis=0)
            ),
            "fc3 supercell frac_coords differ",
        ),
        # matscikit assigns standard masses only, so other masses would change kappa
        (
            lambda fs: fs["unitcell"].update(masses=fs["unitcell"]["masses"] + 1),
            "differs from the standard",
        ),
    ],
)
def test_load_phonons_rejects_inconsistent_force_sets(
    tmp_path: Path,
    cu_calc: ltc.KappaCalculation,
    mutate: Callable[[dict[str, Any]], None],
    match: str,
) -> None:
    """Forces are never fitted to displacements other than the ones they came from."""
    force_set = cu_calc.phonons.force_set()
    mutate(force_set)
    path = write_merged(tmp_path, force_set, "mp")
    with pytest.raises(ValueError, match=match):
        force_sets.load_phonons(path, "mp")


def test_force_set_files_reject_missing_materials_and_other_schemas(
    tmp_path: Path, cu_calc: ltc.KappaCalculation
) -> None:
    """Unknown material IDs and files of another schema version are rejected."""
    merged_path = write_merged(tmp_path, cu_calc.phonons.force_set(), "mp-1")
    with pytest.raises(KeyError, match="mp-2"):
        force_sets.load_force_set(merged_path, "mp-2")
    with h5py.File(merged_path, "a") as h5_file:
        h5_file.attrs["schema_version"] = 2
    with pytest.raises(ValueError, match="force-set schema 2"):
        force_sets.material_ids(merged_path)


def write_schema_2(path: str, phonons: AnharmonicPhonons) -> None:
    """Write phonons' force set as one material in the removed schema-2 layout."""
    unit = phonons.unitcell
    with h5py.File(path, "w") as h5_file:
        h5_file.attrs["schema_version"] = 2
        group = h5_file.create_group("mp-1")
        cell = group.create_group("unitcell")
        for key, value in [
            ("cell", unit.lattice),
            ("scaled_positions", unit.frac_coords),
            ("numbers", unit.numbers),
            ("masses", unit.masses),
        ]:
            cell.create_dataset(key, data=np.asarray(value))
        for key, supercell in [
            ("supercell_matrix", phonons.supercell),
            ("phonon_supercell_matrix", phonons.phonon_supercell),
        ]:
            matrix = np.asarray(supercell.lattice) @ np.linalg.inv(unit.lattice)
            group.create_dataset(key, data=np.rint(matrix))
        group.create_dataset("primitive_matrix", data=phonons.primitive_matrix)
        group.attrs.update(
            symprec=1e-5,
            displacement_distance=0.03,
            plus_minus="null",
            displacement_source="matscikit 0.1",
        )
        for order, supercell, displacements, forces in [
            (
                "fc2",
                phonons.phonon_supercell,
                phonons.phonon_displacements,
                phonons.fc2_forces,
            ),
            ("fc3", phonons.supercell, phonons.displacements, phonons.fc3_forces),
        ]:
            geometry = {
                "lattice": np.asarray(supercell.lattice).tolist(),
                "frac_coords": np.asarray(supercell.frac_coords).tolist(),
            }
            group.attrs[f"{order}_dataset"] = json.dumps(
                {"supercell": geometry, "displacements": displacements}
            )
            group.create_dataset(f"{order}_forces", data=forces.astype(np.float32))


@pytest.mark.parametrize(
    "fc2_supercell",
    [np.diag([2, 2, 2]), np.array([[2, 1, 0], [0, 2, 0], [0, 0, 2]])],
    ids=["diagonal", "sheared"],
)
def test_migrate_force_sets_v2_to_v3(
    tmp_path: Path, cu_calc: ltc.KappaCalculation, fc2_supercell: np.ndarray
) -> None:
    """Migrated schema-2 files keep their forces bit for bit and refit exactly what the
    schema-2 loader fitted (its datasets set verbatim, float32 forces as float64), also
    for a non-symmetric supercell matrix, whose stored form is the transpose.
    """
    from scripts.migrate_force_sets_v2_to_v3 import migrate

    calc = ltc.init_phonons(
        cu_calc.phonons.unitcell,
        fc2_supercell=fc2_supercell,
        fc3_supercell=np.diag([3, 3, 3]),
        q_point_mesh=(4, 4, 4),
        displacement_distance=0.03,
    )
    StandardKappaAdapter().calculate_fc2(calc, EMT(), KappaSettings())
    StandardKappaAdapter().calculate_fc3(calc, EMT(), KappaSettings())
    phonons = calc.phonons
    old_path, new_path = f"{tmp_path}/v2.h5", f"{tmp_path}/v3.h5"
    write_schema_2(old_path, phonons)
    migrate(old_path, new_path)
    with h5py.File(old_path) as old, h5py.File(new_path) as new:
        assert new.attrs["schema_version"] == force_sets.FORCE_SET_SCHEMA_VERSION
        assert new["mp-1"].attrs["displacement_source"] == "matscikit 0.1"
        fc3 = new["mp-1/fc3"]
        np.testing.assert_array_equal(
            np.concatenate([fc3["forces"][()], fc3["pair_forces"][()]]),
            old["mp-1/fc3_forces"][()],
        )
        np.testing.assert_array_equal(
            new["mp-1/fc2/forces"][()], old["mp-1/fc2_forces"][()]
        )
    migrated = force_sets.load_phonons(new_path, "mp-1")
    schema_2_fit = AnharmonicPhonons(
        phonons.unitcell,
        np.diag([3, 3, 3]),
        phonon_supercell_matrix=fc2_supercell,
        primitive_matrix=phonons.primitive_matrix,
    )
    schema_2_fit.set_phonon_displacements(phonons.phonon_displacements, distance=0.03)
    schema_2_fit.set_displacements(phonons.displacements, distance=0.03)
    for order in ("fc2", "fc3"):
        forces = getattr(phonons, f"{order}_forces").astype(np.float32)
        getattr(schema_2_fit, f"produce_{order}")(forces.astype(float))
    np.testing.assert_array_equal(migrated.fc2, schema_2_fit.fc2)
    np.testing.assert_array_equal(migrated.fc3_dense, schema_2_fit.fc3_dense)
    assert migrated.symmetrization == schema_2_fit.symmetrization
    with h5py.File(new_path) as new:  # the measured values are stored, not a marker
        stored = tuple(
            new[f"mp-1/{order}/symmetrization"][()] for order in ("fc2", "fc3")
        )
    assert stored == schema_2_fit.symmetrization
    assert min(stored) > 0
    with pytest.raises(ValueError, match="schema 3, expected 2"):
        migrate(new_path, f"{tmp_path}/again.h5")


def test_migrate_force_sets_v2_to_v3_rejects_unordered_fc3(
    tmp_path: Path, cu_calc: ltc.KappaCalculation
) -> None:
    """FC3 forces split positionally, so a pair config before a single raises."""
    from scripts.migrate_force_sets_v2_to_v3 import migrate

    old_path = f"{tmp_path}/v2.h5"
    write_schema_2(old_path, cu_calc.phonons)
    with h5py.File(old_path, "r+") as old:
        dataset = json.loads(old["mp-1"].attrs["fc3_dataset"])
        configs = dataset["displacements"]
        configs.insert(0, configs.pop())  # move the last pair to the front
        old["mp-1"].attrs["fc3_dataset"] = json.dumps(dataset)
    with pytest.raises(ValueError, match=r"single then .* pair configs"):
        migrate(old_path, f"{tmp_path}/v3.h5")


def test_kappa_runs_save_forces_by_default() -> None:
    """New kappa runs record raw force sets unless a model opts out."""
    assert KappaSettings().save_forces is True
