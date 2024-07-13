import os
import sys
from typing import Any, Literal

import numpy as np
import pandas as pd
import torch
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from pymatgen.core import Structure
#from pymatgen.core.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm import tqdm

from matbench_discovery import timestamp
from matbench_discovery.data import DATA_FILES, as_dict_handler
from matbench_discovery.enums import Key, Task

from sevenn.sevennet_calculator import SevenNetCalculator


__author__ = "Yutack Park"
__date__ = "2024-06-25"


#########################  EDITABLE  ###########################
pot_name = "sevennet"
sevennet_root = None  # root to SevenNet repo
sevennet_checkpoint = f"{sevennet_root}/pretrained_potentials/SevenNet_0__11July2024/checkpoint_sevennet_0.pth"
assert os.path.isfile(sevennet_checkpoint)
task_type = Task.IS2RE
ase_optimizer = "FIRE"
device = "cuda" if torch.cuda.is_available() else "cpu"
ase_filter: Literal["frechet", "exp"] = "frechet"

max_steps = 500
force_max = 0.05  # Run until the forces are smaller than this in eV/A

slurm_array_task_count = 32
#########################  EDITABLE  ###########################
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))

out_dir = "./results"
os.makedirs(out_dir, exist_ok=True)
out_path = f"{out_dir}/{pot_name}-{slurm_array_task_id:>03}.json.gz"

data_path = {
    Task.IS2RE: DATA_FILES.wbm_initial_structures,
}[task_type]
print(f"\nJob started running {timestamp}, eval {pot_name}", flush=True)
print(f"{data_path=}", flush=True)

e_pred_col = "sevennet_energy"

# Init ASE SevenNet Calculator from checkpoint
sevennet_calc = SevenNetCalculator(sevennet_checkpoint)

print(f"Read data from {data_path}")
df_in = pd.read_json(data_path).set_index(Key.mat_id)
df_in = df_in.sample(frac=1, random_state=7)  # shuffle data for equal runtime
if slurm_array_task_count > 1:
    df_in = np.array_split(df_in, slurm_array_task_count)[slurm_array_task_id - 1]

relax_results: dict[str, dict[str, Any]] = {}
input_col = {Task.IS2RE: Key.init_struct}[task_type]

structs = df_in[input_col].map(Structure.from_dict).to_dict()
filter_cls = {"frechet": FrechetCellFilter, "exp": ExpCellFilter}[ase_filter]
optim_cls = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]

for material_id in tqdm(structs, desc="Relaxing"):
    if material_id in relax_results:
        continue
    try:
        atoms = structs[material_id].to_ase_atoms()
        atoms.calc = sevennet_calc
        if max_steps > 0:
            atoms = filter_cls(atoms)
            optimizer = optim_cls(atoms, logfile="/dev/null")
            optimizer.run(fmax=force_max, steps=max_steps)
        energy = atoms.get_potential_energy()  # relaxed energy
        relaxed = AseAtomsAdaptor.get_structure(
            getattr(atoms, "atoms", atoms)  # atoms might be wrapped in ase filter
        )
        relax_results[material_id] = {"structure": relaxed, "energy": energy}

        coords, lattices = (locals().get(key, []) for key in ("coords", "lattices"))
    except Exception as exc:
        print(f"Failed to relax {material_id}: {exc!r}")
        continue

df_out = pd.DataFrame(relax_results).T.add_prefix("sevennet_")
df_out.index.name = Key.mat_id

df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)

