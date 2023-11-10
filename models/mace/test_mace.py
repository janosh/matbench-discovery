# %%
from __future__ import annotations

import os
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
import torch
import wandb
from ase.constraints import ExpCellFilter
from ase.optimize import FIRE, LBFGS
from mace.calculators.mace import MACECalculator
from mace.tools import count_parameters
from pymatgen.core import Structure
from pymatgen.core.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm import tqdm

from matbench_discovery import ROOT, id_col, timestamp, today
from matbench_discovery.data import DATA_FILES, as_dict_handler, df_wbm
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2023-03-01"


# %%
task_type = "IS2RE"  # "RS2RE"
module_dir = os.path.dirname(__file__)
# set large job array size for smaller data splits and faster testing/debugging
slurm_array_task_count = 100
ase_optimizer = "FIRE"
job_name = f"mace-wbm-{task_type}-{ase_optimizer}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")
device = "cuda" if torch.cuda.is_available() else "cpu"
relax_cell = True
model_name = [
    # MACE trained on M3GNet training set by original MACE authors
    "2023-07-14-mace-ilyes-MPF-2021-2-8-big-128-6",
    # MACE trained by Yuan Chiang on CHGNet training set
    "2023-08-14-mace-yuan-mptrj-04",
    "2023-09-03-mace-yuan-mptrj-slower-14-lr-13_run-3",
    "2023-10-29-mace-pbenner-mptrj-no-conditional-loss",
][-1]

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    account="matgen",
    time="4:55:0",
    array=f"1-{slurm_array_task_count}",
    # slurm_flags="--qos shared --constraint gpu --gpus 1",
    slurm_flags="--qos shared --constraint cpu --mem 16G",
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
out_path = f"{out_dir}/mace-preds-{slurm_array_task_id:>03}.json.gz"

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
mace_calc = MACECalculator(checkpoint, device=device)

df_in: pd.DataFrame = np.array_split(
    pd.read_json(data_path).set_index(id_col), slurm_array_task_count
)[slurm_array_task_id - 1]


# %%
run_params = dict(
    data_path=data_path,
    versions={dep: version(dep) for dep in ("mace", "numpy", "torch")},
    checkpoint=checkpoint,
    task_type=task_type,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    slurm_vars=slurm_vars,
    max_steps=max_steps,
    relax_cell=relax_cell,
    force_max=force_max,
    ase_optimizer=ase_optimizer,
    device=device,
    trainable_params=count_parameters(mace_calc.models[0]),
)

run_name = f"{job_name}-{slurm_array_task_id}"
wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %%
relax_results: dict[str, dict[str, Any]] = {}
input_col = {"IS2RE": "initial_structure", "RS2RE": "relaxed_structure"}[task_type]

if task_type == "RS2RE":
    df_in[input_col] = [x["structure"] for x in df_in.computed_structure_entry]

structs = df_in[input_col].map(Structure.from_dict).to_dict()

for material_id in tqdm(structs, desc="Relaxing", disable=None):
    if material_id in relax_results:
        continue
    try:
        atoms = structs[material_id].to_ase_atoms()
        atoms.calc = mace_calc
        if relax_cell:
            atoms = ExpCellFilter(atoms)
            optim_cls = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]
            optimizer = optim_cls(atoms, logfile="/dev/null")

            coords, lattices = [], []
            # attach observer functions to the optimizer
            optimizer.attach(lambda: coords.append(atoms.get_positions()))  # noqa: B023
            optimizer.attach(lambda: lattices.append(atoms.get_cell()))  # noqa: B023

            optimizer.run(fmax=force_max, steps=max_steps)
            mace_traj = Trajectory(
                species=structs[material_id].species,
                coords=coords,
                lattice=lattices,
                constant_lattice=False,
            )
        else:
            mace_traj = None
        mace_energy = atoms.get_potential_energy()
        mace_struct = AseAtomsAdaptor.get_structure(
            atoms.atoms if relax_cell else atoms
        )

        relax_results[material_id] = {
            "structure": mace_struct,
            "energy": mace_energy,
            "trajectory": mace_traj,
        }
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
        ["uncorrected_energy", e_pred_col, "formula"]
    ].reset_index()
)

title = f"MACE {task_type} ({len(df_out):,})"
wandb_scatter(table, fields=dict(x="uncorrected_energy", y=e_pred_col), title=title)

wandb.log_artifact(out_path, type=f"mace-wbm-{task_type}")
