import copy
import json
import os

import numpy as np
from ase.io import read, write
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm import tqdm

from matbench_discovery import FIGSHARE

__author__ = "Yuan Chiang"
__date__ = "2023-08-10"

module_dir = os.path.dirname(__file__)
in_dir = f"{FIGSHARE}/mptrj_2022.9_full.json"
out_dir = f"{module_dir}/mptrj-2022.9"
os.makedirs(out_dir, exist_ok=True)
combined = []


with open(in_dir) as json_file:
    json_data = json.load(json_file)

for material_id in tqdm(json_data):
    for trajectory_id in json_data[material_id]:
        out_xyz = f"{out_dir}/{material_id}.extxyz"

        if os.path.isfile(out_xyz):
            traj = read(out_xyz, index=":", format="extxyz")
            combined.append(traj)
            continue

        block = copy.deepcopy(json_data[material_id][trajectory_id])
        try:
            structure = Structure.from_dict(block.pop("structure"))

            forces = block.pop("force", None)
            magmoms = block.pop("magmom", None)
            stress = block.pop("stress", None)
            # bandgap = block.pop('bandgap', None)

            uncorrected_total_energy = block.pop("uncorrected_total_energy", None)
            mp_id = block.get("mp_id", None)

            atoms = AseAtomsAdaptor.get_atoms(structure=structure)

            if forces:
                atoms.arrays["forces"] = np.array(forces)
            if magmoms:
                atoms.arrays["magmoms"] = np.array(magmoms)
            # if bandgap: will go into atoms.info
            #     atoms.set_tensor('bandgap', bandgap)
            if stress:
                atoms.info["stress"] = np.array(stress)
            if uncorrected_total_energy:
                atoms.info["energy"] = uncorrected_total_energy

            for key, value in block.items():
                atoms.info[key] = value

            assert mp_id == material_id

            write(out_xyz, atoms, append=True, format="extxyz")

            traj = read(out_xyz, index=":", format="extxyz")
            combined.append(traj)

        except Exception as exc:
            print(exc, f"skipping {material_id}, {trajectory_id}")

write("mptrj-2022.9.xyz", combined, format="extxyz")
