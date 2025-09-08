# /// script
# dependencies = [
#   "torch>=2.1.2",
#   "torch-geometric>=2.6.1",
#   "numpy>=1.26.4",
#   "ase>=3.25.0",
#   "braceexpand>=0.1.7",
#   "e3nn>=0.5.6",
#   "pymatviz>=0.16.0",
#   "pyyaml>=6.0.1",
#   "torch-scatter>=2.1.2",
#   "scikit-learn>=1.7.0",
#   "pymatgen>=2025.6.14",
#   "wandb>=0.20.1",
#   "torch-ema>=0.3",
# ]
# ///

import os
from typing import Any, Literal

import pandas as pd
from ase.filters import ExpCellFilter, Filter, FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from hienet.hienet_calculator import HIENetCalculator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import WBM_DIR
from matbench_discovery.data import DataFiles, as_dict_handler, ase_atoms_from_zip
from matbench_discovery.enums import Task

__author__ = "Yutack Park"
__date__ = "2024-06-25"


import argparse

parser = argparse.ArgumentParser(
    description="Run calculations with GPU and data bounds selection"
)
parser.add_argument("--gpu", type=str, default="0", help="GPU ID to use")
parser.add_argument("--left", type=int, default=0, help="Start index for atoms list")
parser.add_argument("--right", type=int, default=None, help="End index for atoms list")
args = parser.parse_args()


# %% this config is editable
left = args.left
right = args.right if args.right != -1 else None


name = "HIENet-V3"
model_name = f"./{name}.pth"
zip_filename = f"{WBM_DIR}/2024-08-04-wbm-initial-atoms.extxyz.zip"
os.makedirs(out_dir := "./results", exist_ok=True)


task_type = Task.IS2RE
ase_optimizer = "FIRE"
os.environ["CUDA_VISIBLE_DEVICES"] = f"{args.gpu}"
device = "cpu"  # "cuda" if torch.cuda.is_available() else "cpu"
ase_filter: Literal["frechet", "exp"] = "frechet"

max_steps = 500
force_max = 0.05  # Run until the forces are smaller than this in eV/A
out_path = f"{out_dir}/{name}_{force_max}_{max_steps}_{left}_{right}.json.gz"

data_path = {Task.IS2RE: DataFiles.wbm_initial_atoms.path}[task_type]

# Initialize ASE Calculator from checkpoint
print(model_name)
hienet_calc = HIENetCalculator(model=model_name)


print(f"Read data from {data_path}")
atoms_list = ase_atoms_from_zip(zip_filename)
atoms_list = atoms_list[left:right]

relax_results: dict[str, dict[str, Any]] = {}

filter_cls: type[Filter] = {
    "frechet": FrechetCellFilter,
    "exp": ExpCellFilter,
}[ase_filter]
optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]


for atoms in tqdm(atoms_list, desc="Relaxing"):
    mat_id = atoms.info[Key.mat_id]
    if mat_id in relax_results:
        continue
    try:
        atoms.calc = hienet_calc

        if max_steps > 0:
            atoms = filter_cls(atoms)
            optimizer = optim_cls(atoms, logfile="/dev/null")
            optimizer.run(fmax=force_max, steps=max_steps)
        energy = atoms.get_potential_energy()  # relaxed energy
        # if max_steps > 0, atoms is wrapped by filter_cls, so extract with getattr
        relaxed_struct = AseAtomsAdaptor.get_structure(getattr(atoms, "atoms", atoms))
        relax_results[mat_id] = {"structure": relaxed_struct, "energy": energy}
    except Exception as exc:
        print(f"Failed to relax {mat_id}: {exc!r}")
        continue

df_out = pd.DataFrame(relax_results).T.add_prefix("hienet_")
df_out.index.name = Key.mat_id

df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)
