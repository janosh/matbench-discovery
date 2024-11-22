# %%
import os
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

from collections.abc import Callable
from copy import deepcopy
from importlib.metadata import version
from typing import Any, Literal

import numpy as np
import pandas as pd
from ase import Atoms
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from pymatgen.core.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm
import json
from matbench_discovery import ROOT, timestamp, today
from matbench_discovery.data import (
    DataFiles,
    as_dict_handler,
    ase_atoms_from_zip,
    df_wbm,
)
from matbench_discovery.enums import MbdKey, Task
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit


from tensorpotential.calculator import grace_fm

__author__ = "Yury Lysogorskiy"
__date__ = "2024-11-21"

# %%
model_name = "MP_GRACE_2L_r6_11Nov2024"

smoke_test = False
task_type = Task.IS2RE
module_dir = os.path.dirname(__file__)
# set large job array size for smaller data splits and faster testing/debugging
slurm_array_task_count = int(os.getenv("SLURM_ARRAY_TASK_COUNT", "1")) # will be set to the number of tasks in the job array.
ase_optimizer = "FIRE"
job_name = f"{model_name}-wbm-{task_type}-{ase_optimizer}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")
device = "gpu"
# whether to record intermediate structures into pymatgen Trajectory
record_traj = False  # has no effect if relax_cell is False

ase_filter: Literal["frechet", "exp"] = "frechet"
os.makedirs(out_dir, exist_ok=True)

# slurm_vars = slurm_submit(
#     job_name=job_name,
#     out_dir=out_dir,
#     array=f"1-{slurm_array_task_count}",
#     # slurm_flags="--qos shared --constraint gpu --gpus 1",
#     slurm_flags="--ntasks=1 --cpus-per-task=1 --partition high-priority",
# )

slurm_vars = {k:v for k,v in os.environ.items() if k.startswith("SLURM_")}

# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0")) #  will be set to the job array index value.
slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug") #  will be set to the first job ID of the array.

out_path = f"{out_dir}/{slurm_array_job_id}-{slurm_array_task_id:>03}.json.gz"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exciting early")

print(f"{slurm_array_task_id=}")
print(f"{slurm_array_job_id=}")
print(f"{slurm_array_task_count=}")
print(f"{slurm_vars=}")
print(f"{out_dir=}")

# %%
data_path = {
    Task.RS2RE: DataFiles.wbm_relaxed_atoms.path,
    Task.IS2RE: DataFiles.wbm_initial_atoms.path,
}[task_type]
print(f"\nJob {job_name} started {timestamp}")
e_pred_col = "grace2l_r6_energy"
max_steps = 500
force_max = 0.05  # Run until the forces are smaller than this in eV/A
checkpoint = "11Nov2024"
dtype = "float64"
calc = grace_fm(model=model_name, pad_neighbors_fraction = 0.05,
                     pad_atoms_number = 2,
                     min_dist=0.5)

print(f"Read data from {data_path}")
atoms_list: list[Atoms] = ase_atoms_from_zip(data_path)
atoms_list=np.array(atoms_list, dtype='object')

if slurm_array_job_id == "debug":
    if smoke_test:
        atoms_list = atoms_list[:128]
    else:
        pass
elif slurm_array_task_count > 1:
    atoms_list = np.array_split(atoms_list, slurm_array_task_count)[slurm_array_task_id - 1]


# %%
run_params = {
    "data_path": data_path,
    "versions": {dep: version(dep) for dep in ("tensorpotential", "numpy", "tensorflow", "ase")},
    "checkpoint": checkpoint,
    Key.task_type: task_type,
    "n_structures": len(atoms_list),
    "slurm_vars": slurm_vars,
    "max_steps": max_steps,
    "record_traj": record_traj,
    "force_max": force_max,
    "ase_optimizer": ase_optimizer,
    "device": device,
    # Key.model_params: count_parameters(calc.models[0]),
    "model_name": model_name,
    "dtype": dtype,
    "ase_filter": ase_filter,
}

run_name = f"{job_name}-{slurm_array_task_id}"

with open(os.path.join(out_dir,f"run_data_{slurm_array_task_id}-{slurm_array_task_count}.json"),"w") as f:
    json.dump(run_params, f)

# wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %% time
relax_results: dict[str, dict[str, Any]] = {}
filter_cls: Callable[[Atoms], Atoms] = {
    "frechet": FrechetCellFilter,
    "exp": ExpCellFilter,
}[ase_filter]
optim_cls: Optimizer = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]
atoms_list = sorted(atoms_list, key=lambda at: len(at))
# print(atoms_list)
for atoms in tqdm(deepcopy(atoms_list), desc="Relaxing", mininterval=5):
    mat_id = atoms.info[Key.mat_id]
    if mat_id in relax_results:
        continue
    try:
        atoms.calc = calc
        if max_steps > 0:
            filtered_atoms = filter_cls(atoms)
            optimizer = optim_cls(filtered_atoms, logfile="/dev/null")

            if record_traj:
                coords, lattices, energies = [], [], []
                # attach observer functions to the optimizer
                optimizer.attach(lambda: coords.append(atoms.get_positions()))  # noqa: B023
                optimizer.attach(lambda: lattices.append(atoms.get_cell()))  # noqa: B023
                optimizer.attach(lambda: energies.append(atoms.get_potential_energy()))  # noqa: B023

            optimizer.run(fmax=force_max, steps=max_steps)
        energy = atoms.get_potential_energy()  # relaxed energy
        # if max_steps > 0, atoms is wrapped by filter_cls, so extract with getattr
        relaxed_struct = AseAtomsAdaptor.get_structure(atoms)
        relax_results[mat_id] = {"structure": relaxed_struct, "energy": energy}

        coords = locals().get("coords", [])
        lattices = locals().get("lattices", [])
        energies = locals().get("energies", [])
        if record_traj and coords and lattices and energies:
            traj = Trajectory(
                species=atoms.get_chemical_symbols(),
                coords=coords,
                lattice=lattices,
                constant_lattice=False,
                frame_properties=[{"energy": energy} for energy in energies],
            )
            relax_results[mat_id]["trajectory"] = traj
    except Exception as exc:
        print(f"Failed to relax {mat_id}: {exc!r}")
        continue


# %%
df_out = pd.DataFrame(relax_results).T.add_prefix("grace2l_r6_")
df_out.index.name = Key.mat_id
if not smoke_test:
    df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)


