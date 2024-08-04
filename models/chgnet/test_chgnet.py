"""Get CHGNet formation energy predictions on WBM test set.
To slurm submit this file: python path/to/file.py slurm-submit
Requires git cloning and then pip installing CHGNet from source:
git clone https://github.com/CederGroupHub/chgnet
pip install -e ./chgnet.
"""

# %%
import os
from importlib.metadata import version
from typing import Any, Literal

import numpy as np
import pandas as pd
import torch
import wandb
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import timestamp, today
from matbench_discovery.data import DataFiles, as_dict_handler, df_wbm
from matbench_discovery.enums import MbdKey, Task
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2023-03-01"

task_type = Task.IS2RE
module_dir = os.path.dirname(__file__)
# set large job array size for smaller data splits and faster testing/debugging
slurm_array_task_count = 50
device = "cuda" if torch.cuda.is_available() else "cpu"
chgnet = StructOptimizer(use_device=device)  # load default pre-trained CHGNnet model
job_name = f"chgnet-{chgnet.version}-wbm-{task_type}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    account="matgen",
    time="11:55:0",
    array=f"1-{slurm_array_task_count}",
    slurm_flags="--qos shared --constraint cpu --nodes 1 --mem 10G",
    # slurm_flags="--qos regular --constraint gpu --gpus 1",
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")
out_path = f"{out_dir}/{slurm_array_job_id}-{slurm_array_task_id:>03}.json.gz"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exciting early")


# %%
data_path = {
    Task.RS2RE: DataFiles.wbm_computed_structure_entries.path,
    Task.IS2RE: DataFiles.wbm_initial_structures.path,
}[task_type]
print(f"\nJob {job_name} started {timestamp}")
print(f"{data_path=}")
e_pred_col = "chgnet_energy"
ase_filter: Literal["FrechetCellFilter", "ExpCellFilter"] = "FrechetCellFilter"
max_steps = 500
fmax = 0.05

df_in = pd.read_json(data_path).set_index(Key.mat_id)
if slurm_array_task_count > 1:
    df_in = np.array_split(df_in, slurm_array_task_count)[slurm_array_task_id - 1]

run_params = {
    "data_path": data_path,
    "versions": {dep: version(dep) for dep in ("chgnet", "numpy", "torch")},
    Key.task_type: task_type,
    "df": {"shape": str(df_in.shape), "columns": ", ".join(df_in)},
    "slurm_vars": slurm_vars,
    "max_steps": max_steps,
    "fmax": fmax,
    "device": device,
    Key.model_params: chgnet.n_params,
}

run_name = f"{job_name}-{slurm_array_task_id}"
wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %%
relax_results: dict[str, dict[str, Any]] = {}
input_col = {Task.IS2RE: Key.init_struct, Task.RS2RE: Key.final_struct}[task_type]

if task_type == Task.RS2RE:
    df_in[input_col] = [cse["structure"] for cse in df_in[Key.cse]]

structures = df_in[input_col].map(Structure.from_dict).to_dict()

for material_id in tqdm(structures, desc="Relaxing"):
    if material_id in relax_results:
        continue
    try:
        relax_result = chgnet.relax(
            structures[material_id],
            verbose=False,
            steps=max_steps,
            fmax=fmax,
            relax_cell=max_steps > 0,
            ase_filter=ase_filter,
        )
        relax_results[material_id] = {
            e_pred_col: relax_result["trajectory"].energies[-1]
        }
        if max_steps > 0:
            relax_struct = relax_result["final_structure"]
            relax_results[material_id]["chgnet_structure"] = relax_struct
            # traj = relax_result["trajectory"]
            # relax_results[material_id]["chgnet_trajectory"] = traj.__dict__
    except Exception as exc:
        print(f"Failed to relax {material_id}: {exc!r}")


# %%
df_out = pd.DataFrame(relax_results).T
df_out.index.name = Key.mat_id

if max_steps == 0:
    df_out.add_suffix("_no_relax").to_csv(out_path.replace(".json.gz", ".csv.gz"))
else:
    df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)


# %%
df_wbm[e_pred_col] = df_out[e_pred_col]
table = wandb.Table(
    dataframe=df_wbm[[MbdKey.dft_energy, e_pred_col, Key.formula]]
    .reset_index()
    .dropna()
)

title = f"CHGNet {task_type} ({len(df_out):,})"
wandb_scatter(table, fields=dict(x=MbdKey.dft_energy, y=e_pred_col), title=title)

wandb.log_artifact(out_path, type=f"chgnet-wbm-{task_type}")
