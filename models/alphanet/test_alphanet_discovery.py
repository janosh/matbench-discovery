import os
from collections.abc import Callable
from typing import Any, Literal

import pandas as pd
import torch
from alphanet.config import All_Config
from alphanet.infer.calc import AlphaNetCalculator
from ase import Atoms
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import timestamp
from matbench_discovery.data import as_dict_handler
from matbench_discovery.enums import Task

smoke_test = False
model_name = "alphanet"
config = All_Config().from_json("./mp/mp.json")

task_type = Task.IS2RE
job_name = f"{model_name}-wbm-{task_type}"
ase_optimizer = "FIRE"
device = "cuda" if torch.cuda.is_available() else "cpu"
ase_filter: Literal["frechet", "exp"] = "frechet"

max_steps = 500
force_max = 0.05  # Run until the forces are smaller than this in eV/A

idx = 1  # we split initial structures into several parts

os.makedirs(out_dir := "./res_relax", exist_ok=True)
out_path = f"{out_dir}/{model_name}-{idx:>03}.json.gz"

data_path = f"./split_relax/part_{idx}.extxyz"
print(f"\nJob {job_name!r} running {timestamp}", flush=True)
print(f"{data_path=}", flush=True)

A_calc = AlphaNetCalculator(
    ckpt_path="./mp/mp_0329.ckpt",
    device="cuda",
    precision="32",
    config=config,
)

print(f"Read data from {data_path}")
atoms_list = read(data_path, index=":", format="extxyz")
relax_results: dict[str, dict[str, Any]] = {}

filter_cls: Callable[[Atoms], Atoms] = {
    "frechet": FrechetCellFilter,
    "exp": ExpCellFilter,
}[ase_filter]
optim_cls: Optimizer = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]

for atoms in tqdm(atoms_list, desc="Relaxing"):
    mat_id = atoms.info[Key.mat_id]
    if mat_id in relax_results:
        continue
    try:
        atoms.calc = A_calc
        if max_steps > 0:
            atoms = filter_cls(atoms)

            optimizer = optim_cls(atoms, logfile="/dev/null")
            optimizer.run(fmax=force_max, steps=max_steps)
        energy = atoms.get_potential_energy()  # relaxed energy
        # if max_steps > 0, atoms is wrapped by filter_cls, so extract with getattr
        relaxed_struct = AseAtomsAdaptor.get_structure(getattr(atoms, "atoms", atoms))
        relax_results[mat_id] = {"structure": relaxed_struct, "energy": energy}
    except Exception:
        print(f"Failed to relax {mat_id}: {exec!r}")

df_out = pd.DataFrame(relax_results).T.add_prefix("alphanet_")
df_out.index.name = Key.mat_id

if not smoke_test:
    df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)
