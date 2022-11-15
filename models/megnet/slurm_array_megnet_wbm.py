# %%
from __future__ import annotations

import os
from datetime import datetime
from importlib.metadata import version

import numpy as np
import pandas as pd
import wandb
from megnet.utils.models import load_model
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.slurm import slurm_submit

"""
To slurm submit this file: python path/to/file.py slurm-submit
Requires Megnet installation: pip install megnet
https://github.com/materialsvirtuallab/megnet
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-11-14"

task_type = "IS2RE"  # "RS2RE"
timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@")[0]
module_dir = os.path.dirname(__file__)
# set large job array size for fast testing/debugging
slurm_array_task_count = 1
slurm_job_id = os.environ.get("SLURM_JOB_ID", "debug")
job_name = f"megnet-wbm-{task_type}-{slurm_job_id}"
out_dir = f"{module_dir}/{today}-{job_name}"

slurm_vars = slurm_submit(
    job_name=job_name,
    log_dir=out_dir,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time=(slurm_max_job_time := "12:0:0"),
    array=f"1-{slurm_array_task_count}",
    # TF_CPP_MIN_LOG_LEVEL=2 means INFO and WARNING logs are not printed
    # https://stackoverflow.com/a/40982782
    pre_cmd="TF_CPP_MIN_LOG_LEVEL=2",
)


# %%
slurm_array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))

print(f"Job started running {timestamp}")

json_out_path = f"{out_dir}/{slurm_array_task_id}.json.gz"
if os.path.isfile(json_out_path):
    raise SystemExit(f"{json_out_path = } already exists, exciting early")


# %%
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
print(f"Loading from {data_path=}")
df_wbm = pd.read_json(data_path).set_index("material_id")

df_this_job: pd.DataFrame = np.array_split(df_wbm, slurm_array_task_count)[
    slurm_array_task_id - 1
]

megnet_mp_e_form = load_model(model_name := "Eform_MP_2019")

run_params = dict(
    data_path=data_path,
    megnet_version=version("megnet"),
    model_name=model_name,
    task_type=task_type,
    slurm_max_job_time=slurm_max_job_time,
    df=dict(shape=str(df_this_job.shape), columns=", ".join(df_this_job)),
    slurm_vars=slurm_vars,
)
if wandb.run is None:
    wandb.login()

wandb.init(
    project="matbench-discovery",
    name=f"{job_name}-{slurm_array_task_id}",
    config=run_params,
)


# %%
if task_type == "IS2RE":
    from pymatgen.core import Structure

    structures = df_this_job.initial_structure.map(Structure.from_dict)
elif task_type == "RS2RE":
    from pymatgen.entries.computed_entries import ComputedStructureEntry

    df_this_job.cse = df_this_job.cse.map(ComputedStructureEntry.from_dict)
    structures = df_this_job.cse.map(lambda x: x.structure)
else:
    raise ValueError(f"Unknown {task_type = }")

megnet_preds = {}
for material_id, structure in tqdm(structures.items(), disable=None):
    if material_id in megnet_preds:
        continue
    e_form_per_atom = megnet_mp_e_form.predict_structure(structure)[0]
    megnet_preds[material_id] = e_form_per_atom


assert len(megnet_preds) == len(structures) == len(df_this_job)
out_col = "megnet_e_form"
df_this_job[out_col] = pd.Series(megnet_preds)


# %%
df_this_job[out_col].reset_index().to_json(json_out_path)

wandb.log_artifact(json_out_path, type=f"m3gnet-wbm-{task_type}")
