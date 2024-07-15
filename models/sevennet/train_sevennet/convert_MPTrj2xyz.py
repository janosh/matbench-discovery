import json
import random

import ase.io
import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import atomic_numbers
from tqdm import tqdm

KBAR_TO_EVpA3 = 1 / 1602.1766208
filename = "MPtrj_2022.9_full.json"

train_ratio = 0.9
val_ratio = 0.1
test_ratio = 0.0

print(f"Reading {filename} ...", flush=True)
with open(filename) as jfile:
    data = json.load(jfile)


def get_id_train_val_test(
    total_size: int,
    train_ratio: float,
    val_ratio: float,
    test_ratio: float,
    split_seed: int | None = 123,
) -> tuple[list, list, list]:
    """Get train, val, test IDs."""
    if train_ratio + val_ratio + test_ratio > 1.0:
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


def chgnet_to_ase_atoms(datum: dict) -> list[Atoms]:
    atoms_list = []
    for m3gid, dtm in datum.items():
        energy = dtm["uncorrected_total_energy"]
        force = dtm["force"]
        stress = dtm["stress"]
        stress = np.array(
            [
                stress[0][0],
                stress[1][1],
                stress[2][2],
                stress[1][2],
                stress[2][0],
                stress[0][1],
            ]
        )
        stress *= -KBAR_TO_EVpA3  # to eV/Angstrom^3
        # internal stress

        stct = dtm["structure"]
        cell = stct["lattice"]["matrix"]
        sites = stct["sites"]
        species = [atomic_numbers[site["species"][0]["element"]] for site in sites]
        pos = [site["xyz"] for site in sites]

        atom = Atoms(species, pos, cell=cell, pbc=True)
        calc_results = {
            "energy": energy,
            "free_energy": energy,
            "forces": force,
            "stress": stress,
        }
        calculator = SinglePointCalculator(atom, **calc_results)
        atom = calculator.get_atoms()

        mpid = m3gid.split("-")[0] + "-" + m3gid.split("-")[1]
        calc_id = m3gid.split("-")[2]
        ionic_step_id = m3gid.split("-")[3]

        info = {
            "data_from": "MP-CHGNet",
            "material_id": mpid,
            "calc_id": calc_id,
            "ionic_step_id": ionic_step_id,
        }
        for if_key in info_keys:
            info[if_key] = dtm[if_key]
        atom.info = info
        atoms_list.append(atom)
    return atoms_list


dataset = list(data.values())
dataset_train = []
dataset_val = []
dataset_test = []

for idx in tqdm(id_train):
    dataset_train.extend(chgnet_to_ase_atoms(dataset[idx]))
ase.io.write("train.extxyz", dataset_train, "extxyz", append=True)

for idx in tqdm(id_val):
    dataset_val.extend(chgnet_to_ase_atoms(dataset[idx]))
ase.io.write("valid.extxyz", dataset_val, "extxyz", append=True)

for idx in tqdm(id_test):
    dataset_test.extend(chgnet_to_ase_atoms(dataset[idx]))
ase.io.write("test.extxyz", dataset_test, "extxyz", append=True)

print("Done!")
