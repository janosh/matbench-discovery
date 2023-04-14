"""Get CHGNet formation energy predictions on WBM test set.
To slurm submit this file: python path/to/file.py slurm-submit
Requires git cloning and then pip installing CHGNet from source:
git clone https://github.com/CederGroupHub/chgnet
pip install -e ./chgnet.
"""


# %%
from __future__ import annotations

import os
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
import wandb
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import DEBUG, timestamp, today
from matbench_discovery.data import DATA_FILES, as_dict_handler, df_wbm
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2023-03-01"

task_type = "IS2RE"  # "RS2RE"
module_dir = os.path.dirname(__file__)
# set large job array size for smaller data splits and faster testing/debugging
slurm_array_task_count = 100
job_name = f"chgnet-wbm-{task_type}{'-debug' if DEBUG else ''}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    partition="ampere",
    account="LEE-SL3-GPU",
    time="6:0:0",
    array=f"1-{slurm_array_task_count}",
    slurm_flags="--nodes 1 --gpus-per-node 1",
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", 0))
out_path = f"{out_dir}/chgnet-preds-{slurm_array_task_id}.json.gz"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")


# %%
data_path = {
    "RS2RE": DATA_FILES.wbm_computed_structure_entries,
    "IS2RE": DATA_FILES.wbm_initial_structures,
}[task_type]
print(f"\nJob started running {timestamp}")
print(f"{data_path=}")
e_pred_col = "chgnet_energy"
max_steps = 2000

df_in: pd.DataFrame = np.array_split(
    pd.read_json(data_path).set_index("material_id"), slurm_array_task_count
)[slurm_array_task_id - 1]

run_params = dict(
    data_path=data_path,
    **{f"{dep}_version": version(dep) for dep in ("chgnet", "numpy", "torch")},
    task_type=task_type,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    slurm_vars=slurm_vars,
    max_steps=max_steps,
)

run_name = f"{job_name}-{slurm_array_task_id}"
wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %%
chgnet = StructOptimizer()  # load default pre-trained CHGNnet model
relax_results: dict[str, dict[str, Any]] = {}
input_col = {"IS2RE": "initial_structure", "RS2RE": "relaxed_structure"}[task_type]

if task_type == "RS2RE":
    df_in[input_col] = [x["structure"] for x in df_in.computed_structure_entry]

structures = df_in[input_col].map(Structure.from_dict).to_dict()

for material_id in tqdm(structures, disable=None):
    if material_id in relax_results:
        continue
    try:
        relax_result = chgnet.relax(
            structures[material_id], verbose=False, steps=max_steps
        )
    except Exception as error:
        print(f"Failed to relax {material_id}: {error}")
        continue
    relax_results[material_id] = {
        "chgnet_structure": relax_result["final_structure"],
        "chgnet_trajectory": relax_result["trajectory"].__dict__,
        e_pred_col: relax_result["trajectory"].energies[-1],
    }


# %%
df_out = pd.DataFrame(relax_results).T
df_out.index.name = "material_id"

df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)


# %%
df_wbm[e_pred_col] = df_out[e_pred_col]
table = wandb.Table(
    dataframe=df_wbm.dropna()[
        ["uncorrected_energy", e_pred_col, "formula"]
    ].reset_index()
)

title = f"CHGNet {task_type} ({len(df_out):,})"
wandb_scatter(table, fields=dict(x="uncorrected_energy", y=e_pred_col), title=title)

wandb.log_artifact(out_path, type=f"chgnet-wbm-{task_type}")
