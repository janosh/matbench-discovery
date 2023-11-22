import copy
import json
import os
import urllib.request

import ase.io
import ase.units
import numpy as np
from pymatgen.core import Structure
from tqdm import tqdm

__author__ = "Yuan Chiang"
__date__ = "2023-08-10"

module_dir = os.path.dirname(__file__)
mptrj_path = f"{module_dir}/MPtrj_2022.9_full.json"
# MPtrj figshare URL
# https://figshare.com/articles/dataset/23713842
urllib.request.urlretrieve(
    "https://figshare.com/ndownloader/files/41619375", mptrj_path
)

with open(mptrj_path) as file:
    json_data = json.load(file)

out_dir = f"{module_dir}/mptrj-gga-ggapu"
os.makedirs(out_dir, exist_ok=True)
combined = []

for material_id in tqdm(json_data):
    xyz_path = f"{out_dir}/{material_id}.extxyz"

    if os.path.isfile(xyz_path):  # read already converted file
        traj = ase.io.read(xyz_path, index=":", format="extxyz")
        combined.extend(traj)
        continue

    for trajectory_id in json_data[material_id]:
        # copy to since can't modify dict while iterating over it
        block = copy.deepcopy(json_data[material_id][trajectory_id])
        try:
            structure = Structure.from_dict(block.pop("structure"))

            forces = block.pop("force", None)
            magmoms = block.pop("magmom", None)
            stress = block.pop("stress", None)

            uncorrected_total_energy = block.pop("uncorrected_total_energy", None)
            mp_id = block.get("mp_id", None)

            atoms = structure.to_ase_atoms()

            if forces:
                atoms.arrays["forces"] = np.array(forces)
            if magmoms:
                atoms.arrays["magmoms"] = np.array(magmoms)
            if stress:
                # kB to eV/A^3
                atoms.info["stress"] = np.array(stress) * 1e-1 * ase.units.GPa
            if uncorrected_total_energy:
                atoms.info["energy"] = uncorrected_total_energy

            for key, value in block.items():
                atoms.info[key] = value

            assert mp_id == material_id

            ase.io.write(xyz_path, atoms, append=True, format="extxyz")

        except Exception as err:
            print(err)

    traj = ase.io.read(xyz_path, index=":", format="extxyz")
    combined.extend(traj)

ase.io.write("mptrj-gga-ggapu.xyz", combined, format="extxyz", append=True)
