"""Tests for HDF5 FC2/FC3 force-set sidecars."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from phono3py import Phono3py
from phonopy.structure.atoms import PhonopyAtoms

from matbench_discovery.phonons import force_sets
from matbench_discovery.phonons import thermal_conductivity as ltc
from matbench_discovery.phonons.pipeline import KappaSettings


@pytest.fixture(scope="module")
def al_force_set() -> tuple[Phono3py, dict[str, Any]]:
    """EMT FC2/FC3 force sets of primitive Al with fitted reference constants."""
    atoms = bulk("Al", "fcc", a=4.05)
    ph3 = Phono3py(
        PhonopyAtoms(
            symbols=atoms.get_chemical_symbols(),
            cell=atoms.cell[:],
            scaled_positions=atoms.get_scaled_positions(),
        ),
        supercell_matrix=3 * np.eye(3, dtype=int),
        phonon_supercell_matrix=3 * np.eye(3, dtype=int),
        primitive_matrix="auto",
    )
    ph3.generate_displacements(distance=0.03)
    ph3.generate_fc2_displacements(distance=0.03)
    fc2_forces = ltc.calculate_fc2_set(ph3, EMT(), pbar_kwargs={"disable": True})
    force_set = force_sets.force_set_from_phono3py(ph3, fc2_forces)
    fc3_forces = ltc.calculate_fc3_set(ph3, EMT(), pbar_kwargs={"disable": True})
    force_set = force_sets.add_fc3_forces(force_set, ph3, fc3_forces)
    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.produce_fc3(symmetrize_fc3r=True)
    return ph3, force_set


def test_force_set_roundtrip_refits_force_constants(
    tmp_path: Path, al_force_set: tuple[Phono3py, dict[str, Any]]
) -> None:
    """Merged HDF5 force sets rebuild FC2/FC3 to float32 force precision."""
    reference, force_set = al_force_set
    shard_path = f"{tmp_path}/shard.h5"
    force_sets.write_force_set_file(shard_path, force_set)
    merged_path = f"{tmp_path}/2026-11-09-phonons-kappa-103-forces.h5"
    force_sets.merge_force_set_files(merged_path, [("mp-134/x", shard_path)])
    assert force_sets.material_ids(merged_path) == ["mp-134/x"]

    stored = force_sets.load_force_set(merged_path, "mp-134/x")
    # forces live only in the float32 arrays, not again inside the datasets
    assert "forces" not in stored["fc3_dataset"]["first_atoms"][0]
    assert stored["fc3_forces"].dtype == np.float32

    ph3 = force_sets.load_phono3py(merged_path, "mp-134/x")
    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.produce_fc3(symmetrize_fc3r=True)
    # float32 keeps ~7 significant digits of each force; measured refit errors are
    # 5e-8 (FC2) and 1e-6 (FC3) of each tensor's largest element, so 1e-5 leaves 10x
    # margin while still catching a lost or misassigned displacement
    for refit, expected in ((ph3.fc2, reference.fc2), (ph3.fc3, reference.fc3)):
        np.testing.assert_allclose(
            refit, expected, rtol=0, atol=1e-5 * np.abs(expected).max()
        )


def test_force_set_without_fc3_and_validation(
    tmp_path: Path, al_force_set: tuple[Phono3py, dict[str, Any]]
) -> None:
    """FC3-gated materials store FC2 only; malformed sets and files are rejected."""
    _reference, force_set = al_force_set
    fc2_only = {**force_set, "fc3_dataset": None, "fc3_forces": None}
    path = f"{tmp_path}/fc2.h5"
    force_sets.write_force_set_file(path, fc2_only)
    merged_path = f"{tmp_path}/merged.h5"
    force_sets.merge_force_set_files(merged_path, [("mp-1", path)])
    ph3 = force_sets.load_phono3py(merged_path, "mp-1")
    assert ph3.forces is None
    assert ph3.phonon_forces.shape == force_set["fc2_forces"].shape

    with pytest.raises(KeyError, match="mp-2"):
        force_sets.load_force_set(merged_path, "mp-2")
    with pytest.raises(ValueError, match="given together"):
        force_sets.write_force_set_file(path, {**force_set, "fc3_forces": None})
    bad_forces = np.full_like(force_set["fc2_forces"], np.nan)
    with pytest.raises(ValueError, match="finite"):
        force_sets.write_force_set_file(path, {**force_set, "fc2_forces": bad_forces})
    with h5py.File(merged_path, "a") as h5_file:
        h5_file.attrs["schema_version"] = 0
    with pytest.raises(ValueError, match="force-set schema 0"):
        force_sets.material_ids(merged_path)


def test_kappa_runs_save_forces_by_default() -> None:
    """New kappa runs record raw force sets unless a model opts out."""
    assert KappaSettings().save_forces is True
