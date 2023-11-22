import copy
import json
import os
import os.path as osp
import wget

import numpy as np
from tqdm.auto import tqdm

from ase import units
from ase.io import read, write
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

__author__ = "Yuan Chiang"
__date__ = "2023-08-10"

mptrj_path = wget.download("https://figshare.com/ndownloader/files/41619375")
mptrj_path = "./MPtrj_2022.9_full.json"

with open(mptrj_path, "r") as f:
    json_data = json.load(f)

#pretty_json_string = json.dumps(json_data, indent=4, ensure_ascii=False)

#mptrj_pretty_path = osp.join(os.curdir, "mptrj_2022.9_pretty.json")
mptrj_extxyz_prefix = osp.join(os.curdir, "mptrj-gga-ggapu")
os.makedirs(mptrj_extxyz_prefix, exist_ok=True)

#with open(mptrj_pretty_path, "r") as f:
#    json_data = json.load(f)

combined = []

for material_id in tqdm(json_data):
    fout_path = osp.join(mptrj_extxyz_prefix, f"{material_id}.extxyz")

    if osp.exists(fout_path):
        traj = read(fout_path, index=":", format="extxyz")
        combined.extend(traj)
        continue

    for trajectory_id in json_data[material_id]:
        block = copy.deepcopy(json_data[material_id][trajectory_id])
        try:
            structure = Structure.from_dict(block.pop("structure"))

            forces = block.pop("force", None)
            magmoms = block.pop("magmom", None)
            stress = block.pop("stress", None)

            uncorrected_total_energy = block.pop("uncorrected_total_energy", None)
            mp_id = block.get("mp_id", None)

            atoms = AseAtomsAdaptor.get_atoms(structure=structure)

            if forces:
                atoms.arrays["forces"] = np.array(forces)
            if magmoms:
                atoms.arrays["magmoms"] = np.array(magmoms)
            if stress:
                atoms.info["stress"] = (
                    np.array(stress) * 1e-1 * units.GPa
                )  # kB to eV/A^3
            if uncorrected_total_energy:
                atoms.info["energy"] = uncorrected_total_energy

            for key, value in block.items():
                atoms.info[key] = value

            assert mp_id == material_id

            write(fout_path, atoms, append=True, format="extxyz")

        except Exception as err:
            print(err)

    traj = read(fout_path, index=":", format="extxyz")
    combined.extend(traj)

write("mptrj-gga-ggapu.xyz", combined, format="extxyz", append=True)
