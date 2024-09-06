"""Converts energy/force/stress-decorated pymatgen structures to ASE Atoms and
writes them to extxyz files.
"""

# %%
import gzip
import json

import ase.io
import ase.units
import numpy as np
import pandas as pd
from ase import Atoms
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import MP_DIR, WBM_DIR, today
from matbench_discovery.data import DataFiles, ase_atoms_to_zip

__author__ = "Yuan Chiang, Janosh Riebesell"
__date__ = "2023-08-10"

# TODO 2024-08-04 figure out why ase.io.read for Atoms in gzipped extxyz format scales
# superlinearly (0.3 sec for 1_000 Atoms, 39 sec for 10_000, 250k unknown since never
# finished). .zip doesn't have this problem but file size is about 4x larger than .gz.


# %% convert MPtrj pymatgen Structure to ASE Atoms with Structure.properties mapped
# to Atoms.info
with gzip.open(DataFiles.mp_trj.path, mode="rt") as file:
    mptrj_data = json.load(file)

mptrj_atoms_list: list[Atoms] = []
for mat_id, trajectories in tqdm(mptrj_data.items(), desc="MPtrj"):
    for trajectory_id, struct_with_meta in trajectories.items():
        struct = Structure.from_dict(struct_with_meta["structure"])
        struct.properties[Key.mat_id] = mat_id
        struct.properties[Key.frame_id] = trajectory_id
        struct.properties |= struct_with_meta.copy()
        struct.properties.pop("structure")
        atoms = struct.to_ase_atoms()

        if stress := struct.properties.get("stress"):
            struct.properties["stress_kBar"] = stress
            # Convert stress from kBar to eV/A^3 and use ASE sign convention
            atoms.info["stress"] = np.array(stress) * -0.1 * ase.units.GPa
        mptrj_atoms_list.append(atoms)

ase_atoms_to_zip(mptrj_atoms_list, f"{MP_DIR}/{today}-mp-trj.extxyz.zip")


# %% convert WBM initial structures to ASE Atoms (no properties other than material ID
# included in Atoms.info)
df_wbm_init = pd.read_json(DataFiles.wbm_initial_structures.path).set_index(Key.mat_id)

wbm_init_atoms_list: list[Atoms] = []
for mat_id, struct_dict in tqdm(df_wbm_init[Key.init_struct].items(), desc="WBM init"):
    struct = Structure.from_dict(struct_dict)
    struct.properties[Key.mat_id] = mat_id
    atoms = struct.to_ase_atoms()
    wbm_init_atoms_list += [atoms]

ase_atoms_to_zip(wbm_init_atoms_list, f"{WBM_DIR}/{today}-wbm-initial-atoms.extxyz.zip")


# %% convert WBM ComputedStructureEntries to ASE Atoms (material ID and energy included
# in Atoms.info)
df_wbm_cse = pd.read_json(DataFiles.wbm_computed_structure_entries.path).set_index(
    Key.mat_id
)

wbm_cse_atoms_list: list[Atoms] = []
for mat_id, cse_dict in tqdm(df_wbm_cse[Key.cse].items(), desc="WBM CSE"):
    struct = Structure.from_dict(cse_dict[Key.structure])
    struct.properties[Key.mat_id] = mat_id
    struct.properties[Key.energy] = cse_dict[Key.energy]
    atoms = struct.to_ase_atoms()
    wbm_cse_atoms_list += [atoms]

ase_atoms_to_zip(wbm_cse_atoms_list, f"{WBM_DIR}/{today}-wbm-relaxed-atoms.extxyz.zip")
