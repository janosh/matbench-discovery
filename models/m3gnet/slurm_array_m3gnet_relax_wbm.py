# %%
from __future__ import annotations

import os
import warnings
from datetime import datetime
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
import wandb
from m3gnet.models import Relaxer
from tqdm import tqdm

from mb_discovery import ROOT, as_dict_handler
from mb_discovery.slurm import slurm_submit_python

"""
To slurm submit this file, run:

python path/to/file.py slurm-submit

Requires M3GNet installation: pip install m3gnet
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

task_type = "IS2RE"  # "RS2RE"
today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)
# set large job array size for fast testing/debugging
slurm_array_task_count = 100
slurm_mem_per_node = 12000
job_name = f"m3gnet-wbm-relax-{task_type}"
out_dir = f"{module_dir}/{today}-{job_name}"

slurm_submit_python(
    job_name=job_name,
    log_dir=out_dir,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time=(slurm_max_job_time := "3:0:0"),
    array=f"1-{slurm_array_task_count}",
    slurm_flags=("--mem", str(slurm_mem_per_node)),
    # TF_CPP_MIN_LOG_LEVEL=2 means INFO and WARNING logs are not printed
    # https://stackoverflow.com/a/40982782
    pre_cmd="TF_CPP_MIN_LOG_LEVEL=2",
)


# %%
slurm_job_id = os.environ.get("SLURM_JOB_ID", "debug")
slurm_array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))

print(f"Job started running {datetime.now():%Y-%m-%d@%H-%M}")
print(f"{slurm_job_id = }")
print(f"{slurm_array_task_id = }")
print(f"{version('m3gnet') = }")

json_out_path = f"{out_dir}/{slurm_array_task_id}.json.gz"

if os.path.isfile(json_out_path):
    raise SystemExit(f"{json_out_path = } already exists, exciting early")

warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")
warnings.filterwarnings(action="ignore", category=UserWarning, module="tensorflow")


# %%
data_path = f"{ROOT}/data/2022-06-26-wbm-cses-and-initial-structures.json.gz"
print(f"Loading from {data_path=}")
df_wbm = pd.read_json(data_path).set_index("material_id")

df_this_job = np.array_split(df_wbm, slurm_array_task_count)[slurm_array_task_id - 1]

run_params = dict(
    data_path=data_path,
    m3gnet_version=version("m3gnet"),
    slurm_array_task_count=slurm_array_task_count,
    slurm_array_task_id=slurm_array_task_id,
    slurm_job_id=slurm_job_id,
    slurm_max_job_time=slurm_max_job_time,
    slurm_mem_per_node=slurm_mem_per_node,
    task_type=task_type,
)
if wandb.run is None:
    wandb.login()

wandb.init(
    project="m3gnet",
    name=f"{job_name}-{slurm_job_id}-{slurm_array_task_id}",
    config=run_params,
)


# %%
relaxer = Relaxer()  # This loads the default pre-trained M3GNet model
relax_results: dict[str, dict[str, Any]] = {}

if task_type == "IS2RE":
    from pymatgen.core import Structure

    structures = df_this_job.initial_structure.map(Structure.from_dict)
elif task_type == "RS2RE":
    from pymatgen.entries.computed_entries import ComputedStructureEntry

    df_this_job.cse = df_this_job.cse.map(ComputedStructureEntry.from_dict)
    structures = df_this_job.cse.map(lambda x: x.structure)
else:
    raise ValueError(f"Unknown {task_type = }")


for material_id, structure in tqdm(structures.items(), disable=None):
    if material_id in relax_results:
        continue
    relax_result = relaxer.relax(structure)
    relax_dict = {
        "m3gnet_structure": relax_result["final_structure"],
        "m3gnet_trajectory": relax_result["trajectory"].__dict__,
    }

    relax_results[material_id] = relax_dict


# %%
df_output = pd.DataFrame(relax_results).T
df_output.index.name = "material_id"

df_output.reset_index().to_json(json_out_path, default_handler=as_dict_handler)

wandb.log_artifact(json_out_path, type=f"m3gnet-relax-wbm-{task_type}")
