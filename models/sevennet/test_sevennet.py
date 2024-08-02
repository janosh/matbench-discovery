# %%
import os
from typing import Any, Literal

import numpy as np
import pandas as pd
import torch
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from sevenn.sevennet_calculator import SevenNetCalculator
from tqdm import tqdm

from matbench_discovery import timestamp
from matbench_discovery.data import DataFiles, as_dict_handler
from matbench_discovery.enums import Task

__author__ = "Yutack Park"
__date__ = "2024-06-25"


# %% this config is editable
smoke_test = True
model_name = "7net-0"
task_type = Task.IS2RE
ase_optimizer = "FIRE"
device = "cuda" if torch.cuda.is_available() else "cpu"
ase_filter: Literal["frechet", "exp"] = "frechet"

max_steps = 500
force_max = 0.05  # Run until the forces are smaller than this in eV/A

slurm_array_task_count = 32


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))

os.makedirs(out_dir := "./results", exist_ok=True)
out_path = f"{out_dir}/{model_name}-{slurm_array_task_id:>03}.json.gz"

data_path = {Task.IS2RE: DataFiles.wbm_initial_structures.path}[task_type]
print(f"\nJob started running {timestamp}, eval {model_name}", flush=True)
print(f"{data_path=}", flush=True)

# Initialize ASE SevenNet Calculator from checkpoint
seven_net_calc = SevenNetCalculator(model=model_name)


# %%
print(f"Read data from {data_path}")
df_in = pd.read_json(data_path).set_index(Key.mat_id)
if smoke_test:
    df_in = df_in.head(10)
else:
    df_in = df_in.sample(frac=1, random_state=7)  # shuffle data for equal runtime
    if slurm_array_task_count > 1:
        df_in = np.array_split(df_in, slurm_array_task_count)[slurm_array_task_id - 1]

relax_results: dict[str, dict[str, Any]] = {}
input_col = {Task.IS2RE: Key.init_struct}[task_type]

structs = df_in[input_col].map(Structure.from_dict).to_dict()
filter_cls = {"frechet": FrechetCellFilter, "exp": ExpCellFilter}[ase_filter]
optim_cls = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]


# %%
for mat_id in tqdm(structs, desc="Relaxing"):
    if mat_id in relax_results:
        continue
    try:
        atoms = structs[mat_id].to_ase_atoms()
        atoms.calc = seven_net_calc
        if max_steps > 0:
            atoms = filter_cls(atoms)
            optimizer = optim_cls(atoms, logfile="/dev/null")
            optimizer.run(fmax=force_max, steps=max_steps)
        energy = atoms.get_potential_energy()  # relaxed energy
        # atoms might be wrapped in ase filter
        relaxed = AseAtomsAdaptor.get_structure(getattr(atoms, "atoms", atoms))
        relax_results[mat_id] = {"structure": relaxed, "energy": energy}

        coords, lattices = (locals().get(key, []) for key in ("coords", "lattices"))
    except Exception as exc:
        print(f"Failed to relax {mat_id}: {exc!r}")
        continue

df_out = pd.DataFrame(relax_results).T.add_prefix("sevennet_")
df_out.index.name = Key.mat_id


# %%
if not smoke_test:
    df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)
