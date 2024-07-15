#!/usr/bin/env python
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

train_ratio = 0.8
val_ratio = 0.1
test_ratio = 0.1

print(f"Reading {filename} ...", flush=True)
with open(filename) as jfile:
    data = json.load(jfile)


def get_id_train_val_test(
    total_size=1000,
    split_seed=123,
    train_ratio=None,
    val_ratio=0.1,
    test_ratio=0.1,
    n_train=None,
    n_test=None,
    n_val=None,
    keep_data_order=False,
):
    """Get train, val, test IDs."""
    if train_ratio is None and val_ratio is not None and test_ratio is not None:
        if train_ratio is None:
            assert val_ratio + test_ratio < 1
            train_ratio = 1 - val_ratio - test_ratio
            print("Using rest of the dataset except the test and val sets.")
        else:
            assert train_ratio + val_ratio + test_ratio <= 1
    if n_train is None:
        n_train = int(train_ratio * total_size)
    if n_test is None:
        n_test = int(test_ratio * total_size)
    if n_val is None:
        n_val = int(val_ratio * total_size)
    ids = list(np.arange(total_size))
    if not keep_data_order:
        random.seed(split_seed)
        random.shuffle(ids)
    if n_train + n_val + n_test > total_size:
        raise ValueError(
            "Check total number of samples.",
            n_train + n_val + n_test,
            ">",
            total_size,
        )

    id_train = ids[:n_train]
    id_val = ids[-(n_val + n_test) : -n_test]
    id_test = ids[-n_test:]
    return id_train, id_val, id_test


id_train, id_val, id_test = get_id_train_val_test(
    total_size=len(data),
    split_seed=1,
    train_ratio=train_ratio,
    val_ratio=val_ratio,
    test_ratio=test_ratio,
    keep_data_order=False,
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


def chgnet_to_ase_atoms(datum):
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
