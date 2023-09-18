import copy
import json
import os
import os.path as osp

import numpy as np
from ase.io import read, write
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm.auto import tqdm

__author__ = "Yuan Chiang"
__date__ = "2023-08-10"

mptrj_pretty_path = osp.join("../../data/figshare", "mptrj_2022.9_full.json")
mptrj_extxyz_prefix = osp.join(os.curdir, "mptrj-2022.9")

os.makedirs(mptrj_extxyz_prefix, exist_ok=True)

with open(mptrj_pretty_path) as f:
    json_data = json.load(f)

combined = []

for material_id in tqdm(json_data):
    for trajectory_id in json_data[material_id]:
        fout_path = osp.join(mptrj_extxyz_prefix, f"{material_id}.extxyz")

        if osp.exists(fout_path):
            traj = read(fout_path, index=":", format="extxyz")
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

            write(fout_path, atoms, append=True, format="extxyz")

            traj = read(fout_path, index=":", format="extxyz")
            combined.append(traj)

        except Exception as err:
            print(err, f"skipping {material_id}, {trajectory_id}")

write("mptrj-2022.9.xyz", combined, format="extxyz")
