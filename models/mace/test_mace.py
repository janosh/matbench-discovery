# %%
from __future__ import annotations

import os
from importlib.metadata import version
from typing import Any, Literal

import numpy as np
import pandas as pd
import torch
import wandb
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from mace.calculators import mace_mp
from mace.tools import count_parameters
from pymatgen.core import Structure
from pymatgen.core.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm import tqdm

from matbench_discovery import ROOT, formula_col, id_col, timestamp, today
from matbench_discovery.data import DATA_FILES, as_dict_handler, df_wbm
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2023-03-01"


# %%
task_type = "IS2RE"  # "RS2RE"
module_dir = os.path.dirname(__file__)
# set large job array size for smaller data splits and faster testing/debugging
slurm_array_task_count = 50
ase_optimizer = "FIRE"
job_name = f"mace-wbm-{task_type}-{ase_optimizer}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")
device = "cuda" if torch.cuda.is_available() else "cpu"
# whether to record intermediate structures into pymatgen Trajectory
record_traj = False  # has no effect if relax_cell is False
model_name = [
    "2023-10-29-mace-16M-pbenner-mptrj-no-conditional-loss",
    "https://tinyurl.com/y7uhwpje",
][-1]
ase_filter: Literal["frechet", "exp"] = "frechet"

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    account="matgen",
    time="11:55:0",
    array=f"1-{slurm_array_task_count}",
    # slurm_flags="--qos shared --constraint gpu --gpus 1",
    slurm_flags="--qos shared --constraint cpu --mem 32G",
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")
out_path = f"{out_dir}/{slurm_array_job_id}-{slurm_array_task_id:>03}.json.gz"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exciting early")


# %%
data_path = {
    "RS2RE": DATA_FILES.wbm_computed_structure_entries,
    "IS2RE": DATA_FILES.wbm_initial_structures,
}[task_type]
print(f"\nJob started running {timestamp}")
print(f"{data_path=}")
e_pred_col = "mace_energy"
max_steps = 500
force_max = 0.05  # Run until the forces are smaller than this in eV/A
checkpoint = f"{ROOT}/models/mace/checkpoints/{model_name}.model"
dtype = "float64"
mace_calc = mace_mp(model=model_name, device=device, default_dtype=dtype)

df_in = pd.read_json(data_path).set_index(id_col)
if slurm_array_task_count > 1:
    df_in = np.array_split(df_in, slurm_array_task_count)[slurm_array_task_id - 1]


# %%
run_params = dict(
    data_path=data_path,
    versions={dep: version(dep) for dep in ("mace", "numpy", "torch")},
    checkpoint=checkpoint,
    task_type=task_type,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    slurm_vars=slurm_vars,
    max_steps=max_steps,
    record_traj=record_traj,
    force_max=force_max,
    ase_optimizer=ase_optimizer,
    device=device,
    trainable_params=count_parameters(mace_calc.models[0]),
    model_name=model_name,
    dtype=dtype,
    ase_filter=ase_filter,
)

run_name = f"{job_name}-{slurm_array_task_id}"
wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %%
relax_results: dict[str, dict[str, Any]] = {}
input_col = {"IS2RE": "initial_structure", "RS2RE": "relaxed_structure"}[task_type]

if task_type == "RS2RE":
    df_in[input_col] = [x["structure"] for x in df_in.computed_structure_entry]

structs = df_in[input_col].map(Structure.from_dict).to_dict()
filter_cls = {"frechet": FrechetCellFilter, "exp": ExpCellFilter}[ase_filter]

for material_id in tqdm(structs, desc="Relaxing"):
    if material_id in relax_results:
        continue
    try:
        mace_traj = None
        atoms = structs[material_id].to_ase_atoms()
        atoms.calc = mace_calc
        if max_steps > 0:
            atoms = filter_cls(atoms)
            optim_cls = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]
            optimizer = optim_cls(atoms, logfile="/dev/null")

            if record_traj:
                coords, lattices = [], []
                # attach observer functions to the optimizer
                optimizer.attach(lambda: coords.append(atoms.get_positions()))  # noqa: B023
                optimizer.attach(lambda: lattices.append(atoms.get_cell()))  # noqa: B023

            optimizer.run(fmax=force_max, steps=max_steps)
        mace_energy = atoms.get_potential_energy()  # relaxed energy
        mace_struct = AseAtomsAdaptor.get_structure(
            getattr(atoms, "atoms", atoms)  # atoms might be wrapped in ase filter
        )

        relax_results[material_id] = {"structure": mace_struct, "energy": mace_energy}
        if record_traj and len(coords) > 0:
            mace_traj = Trajectory(
                species=structs[material_id].species,
                coords=coords,
                lattice=lattices,
                constant_lattice=False,
            )
            relax_results[material_id]["trajectory"] = mace_traj
    except Exception as exc:
        print(f"Failed to relax {material_id}: {exc!r}")
        continue


# %%
df_out = pd.DataFrame(relax_results).T.add_prefix("mace_")
df_out.index.name = id_col

df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)


# %%
df_wbm[e_pred_col] = df_out[e_pred_col]
table = wandb.Table(
    dataframe=df_wbm.dropna()[
        ["uncorrected_energy", e_pred_col, formula_col]
    ].reset_index()
)

title = f"MACE {task_type} ({len(df_out):,})"
wandb_scatter(table, fields=dict(x="uncorrected_energy", y=e_pred_col), title=title)

wandb.log_artifact(out_path, type=f"mace-wbm-{task_type}")
