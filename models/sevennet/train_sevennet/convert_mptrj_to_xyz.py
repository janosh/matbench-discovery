import json
import random
from typing import Any

import ase.io
import ase.units
import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import atomic_numbers
from ase.stress import full_3x3_to_voigt_6_stress
from pymatviz.enums import Key
from tqdm import tqdm

# Convert stress from kBar to eV/A^3 and use ASE sign convention
kbar_to_evpa3 = -0.1 / ase.units.GPa
filename = "MPtrj_2022.9_full.json"
train_ratio, val_ratio, test_ratio = 0.9, 0.1, 0.0

print(f"Reading {filename} ...", flush=True)
with open(filename) as file:
    data = json.load(file)


def get_id_train_val_test(
    total_size: int,
    train_ratio: float,
    val_ratio: float,
    test_ratio: float,
    split_seed: int = 123,
) -> tuple[list[int], list[int], list[int]]:
    """Get train, val, test IDs."""
    if train_ratio + val_ratio + test_ratio > 1:
        raise ValueError("train_ratio + val_ratio + test_ratio is over 1.0")
    n_train = int(train_ratio * total_size)
    n_val = int(val_ratio * total_size)
    n_test = int(test_ratio * total_size)
    ids = list(np.arange(total_size))

    random.seed(split_seed)
    random.shuffle(ids)

    id_train = ids[:n_train]
    id_val = ids[-(n_val + n_test) : -n_test]
    id_test = ids[-n_test:] if n_test != 0 else []
    return id_train, id_val, id_test


id_train, id_val, id_test = get_id_train_val_test(
    total_size=len(data),
    train_ratio=train_ratio,
    val_ratio=val_ratio,
    test_ratio=test_ratio,
)

info_keys = [
    "uncorrected_total_energy",
    "ef_per_atom",
    "e_per_atom_relaxed",
    "ef_per_atom_relaxed",
    "magmom",
    "bandgap",
    "mp_id",
]


def chgnet_to_ase_atoms(datum: dict[str, dict[str, Any]]) -> list[Atoms]:
    atoms_list = []
    for mat_id, dtm in datum.items():
        energy = dtm["uncorrected_total_energy"]
        force = dtm["force"]
        stress = full_3x3_to_voigt_6_stress(dtm["stress"])  # internal stress
        stress *= kbar_to_evpa3  # to eV/Angstrom^3

        struct = dtm["structure"]
        cell = struct["lattice"]["matrix"]
        sites = struct["sites"]
        species = [atomic_numbers[site["species"][0]["element"]] for site in sites]
        pos = [site["xyz"] for site in sites]

        atoms = Atoms(species, pos, cell=cell, pbc=True)
        calc_results = {
            Key.energy: energy,
            Key.free_energy: energy,
            Key.forces: force,
            Key.stress: stress,
        }
        calculator = SinglePointCalculator(atoms, **calc_results)
        atoms = calculator.get_atoms()

        info = {
            "data_from": "MP-CHGNet",
            "material_id": mat_id.split("-")[0] + "-" + mat_id.split("-")[1],
            "calc_id": mat_id.split("-")[2],
            "ionic_step_id": mat_id.split("-")[3],
        }
        for key in info_keys:
            info[key] = dtm[key]
        atoms.info = info
        atoms_list += [atoms]
    return atoms_list


dataset = list(data.values())
train_set, val_set, test_set = [], [], []

for idx in tqdm(id_train):
    train_set.extend(chgnet_to_ase_atoms(dataset[idx]))
ase.io.write("train.extxyz", train_set, "extxyz", append=True)

for idx in tqdm(id_val):
    val_set.extend(chgnet_to_ase_atoms(dataset[idx]))
ase.io.write("valid.extxyz", val_set, "extxyz", append=True)

for idx in tqdm(id_test):
    test_set.extend(chgnet_to_ase_atoms(dataset[idx]))
ase.io.write("test.extxyz", test_set, "extxyz", append=True)

print("Done!")
