"""
Script for geometry optimization and convex hull calculations, using a PET model.

Adapted from the corresponding GRACE script.
"""

import json
import os
import warnings
from copy import deepcopy
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from metatomic.torch.ase_calculator import MetatomicCalculator, SymmetrizedCalculator
from pymatgen.core.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import timestamp
from matbench_discovery.data import as_dict_handler, ase_atoms_from_zip
from matbench_discovery.enums import DataFiles, Task

warnings.filterwarnings("ignore", category=RuntimeWarning)


model_name = "pet"
smoke_test = False
task_type = Task.IS2RE
module_dir = os.path.dirname(__file__)
# set large job array size for smaller data splits and faster testing/debugging
slurm_array_task_count = int(
    os.getenv("SLURM_ARRAY_TASK_COUNT", "1")
)  # will be set to the number of tasks in the job array.
ase_optimizer = "FIRE"
job_name = "pet"
out_dir = os.getenv("SBATCH_OUTPUT", f"{job_name}")
device = "cuda"
# whether to record intermediate structures into pymatgen Trajectory
record_traj = False  # has no effect if relax_cell is False
os.makedirs(out_dir, exist_ok=True)


# will be set to the job array index value.
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
# will be set to the first job ID of the array.
slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")

out_path = f"{out_dir}/{slurm_array_job_id}-{slurm_array_task_id:>03}.json.gz"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exiting early")

print(f"{slurm_array_task_id=}")
print(f"{slurm_array_job_id=}")
print(f"{slurm_array_task_count=}")
print(f"{out_dir=}")


data_path = {
    Task.RS2RE: DataFiles.wbm_relaxed_atoms.path,
    Task.IS2RE: DataFiles.wbm_initial_atoms.path,
}[task_type]
print(f"\nJob {job_name} started {timestamp}")
e_pred_col = "pet_energy"
max_steps = 500
force_max = 0.02  # Run until the forces are smaller than this in eV/A
checkpoint = ""

# get it with `mtt export https://huggingface.co/lab-cosmo/upet/resolve/main/models/pet-oam-xl-v1.0.0.ckpt`
calc = MetatomicCalculator("pet-oam-xl-v1.0.0.pt", device=device)
calc = SymmetrizedCalculator(calc, batch_size=16, include_inversion=False)

print(f"Read data from {data_path}")
atoms_list = ase_atoms_from_zip(data_path)
atoms_list = np.array(atoms_list, dtype="object")

if slurm_array_job_id == "debug":
    if smoke_test:
        atoms_list = atoms_list[:128]
    else:
        pass
elif slurm_array_task_count > 1:
    atoms_list = np.array_split(atoms_list, slurm_array_task_count)[
        slurm_array_task_id - 1
    ]


run_params = {
    "data_path": data_path,
    "versions": {dep: version(dep) for dep in ("numpy", "ase")},
    "checkpoint": checkpoint,
    Key.task_type: task_type,
    "n_structures": len(atoms_list),
    "max_steps": max_steps,
    "record_traj": record_traj,
    "force_max": force_max,
    "ase_optimizer": ase_optimizer,
    "device": device,
    "model_name": model_name,
    "cell_filter": "FrechetCellFilter",
}

run_name = f"{job_name}-{slurm_array_task_id}"

with open(
    f"{out_dir}/run_data_{slurm_array_task_id}-{slurm_array_task_count}.json", mode="w"
) as file:
    json.dump(run_params, file)


relax_results: dict[str, dict[str, Any]] = {}
optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]
atoms_list = sorted(atoms_list, key=lambda at: len(at))
# print(atoms_list)
for atoms in tqdm(deepcopy(atoms_list), desc="Relaxing", mininterval=5):
    mat_id = atoms.info[Key.mat_id]
    if mat_id in relax_results:
        continue
    try:
        atoms.calc = calc
        if max_steps > 0:
            filtered_atoms = FrechetCellFilter(atoms)
            optimizer = optim_cls(filtered_atoms, logfile="/dev/null")

            if record_traj:
                coords, lattices, energies = [], [], []
                # attach observer functions to the optimizer
                optimizer.attach(lambda: coords.append(atoms.get_positions()))  # noqa: B023
                optimizer.attach(lambda: lattices.append(atoms.get_cell()))  # noqa: B023
                optimizer.attach(
                    lambda: energies.append(atoms.get_potential_energy())  # noqa: B023
                )

            optimizer.run(fmax=force_max, steps=max_steps)
        energy = atoms.get_potential_energy()  # relaxed energy
        # if max_steps > 0, atoms is wrapped by FrechetCellFilter, so need to getattr
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


df_out = pd.DataFrame(relax_results).T.add_prefix("pet_")
df_out.index.name = Key.mat_id
if not smoke_test:
    df_out.reset_index().to_json(
        out_path, default_handler=as_dict_handler, orient="records", lines=True
    )
