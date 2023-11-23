# %%
import json
import os
import re
import urllib.request

import ase.io
import ase.units
import numpy as np
from pymatgen.core import Structure
from pymatviz.io import TqdmDownload
from tqdm import tqdm

from matbench_discovery import DATA_DIR

__author__ = "Yuan Chiang"
__date__ = "2023-08-10"

module_dir = os.path.dirname(__file__)
mp_trj_path = f"{DATA_DIR}/MPtrj_2022.9_full.json"


# %% MPtrj figshare URL https://figshare.com/articles/dataset/23713842
# the download is 11.3 GB and can easily take 1h
mp_trj_url = "https://figshare.com/ndownloader/files/41619375"


with TqdmDownload(desc=mp_trj_url) as pbar:
    urllib.request.urlretrieve(mp_trj_url, mp_trj_path, reporthook=pbar.update_to)


# %%
with open(mp_trj_path) as file:
    json_data = json.load(file)

out_dir = f"{module_dir}/mptrj-gga-ggapu"
os.makedirs(out_dir, exist_ok=True)
combined = []


# %%
for material_id in tqdm(json_data):
    xyz_path = f"{out_dir}/{material_id}.extxyz"

    if os.path.isfile(xyz_path):  # read already converted file
        traj = ase.io.read(xyz_path, index=":", format="extxyz")
        combined.extend(traj)
        continue

    for trajectory_id in json_data[material_id]:
        block = json_data[material_id][trajectory_id]
        try:
            atoms = Structure.from_dict(block["structure"]).to_ase_atoms()

            match = re.match(r"(mp-\d+)-(\d+)-(\d+)", trajectory_id)
            if not match:
                raise ValueError(f"Invalid {trajectory_id=}")

            task_id, calc_id, ionic_step = match.groups()
            atoms.info["task_id"] = task_id
            atoms.info["calc_id"] = int(calc_id)
            atoms.info["ionic_step"] = int(ionic_step)

            mp_id = block.get("mp_id")
            assert mp_id == material_id, f"{mp_id=} != {material_id=}"

            if uncorrected_total_energy := block.get("uncorrected_total_energy"):
                atoms.info["energy"] = uncorrected_total_energy
            if forces := block.get("force"):
                atoms.arrays["forces"] = np.array(forces)
            if magmoms := block.get("magmom"):
                atoms.arrays["magmoms"] = np.array(magmoms)
            if stress := block.get("stress"):
                # kB to eV/A^3 and opposite sign convention
                atoms.info["stress"] = np.array(stress) * -1e-1 * ase.units.GPa

            special_keys = {"uncorrected_total_energy", "force", "magmom", "stress"}
            for key in {*block} - special_keys:
                if val := atoms.get(key):
                    atoms.info[key] = val

            ase.io.write(xyz_path, atoms, append=True, format="extxyz")

        except Exception as err:
            print(err)

    traj = ase.io.read(xyz_path, index=":", format="extxyz")
    combined.extend(traj)


# %%
ase.io.write("mptrj-gga-ggapu.xyz", combined, format="extxyz", append=True)
