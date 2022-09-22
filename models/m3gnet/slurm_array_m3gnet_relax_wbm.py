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

from mb_discovery import ROOT, as_dict_handler

"""
To slurm submit this file, use

```sh
# slurm will not create logdir automatically and fail if missing
mkdir -p models/m3gnet/slurm_logs
sbatch --partition icelake-himem --account LEE-SL3-CPU --array 1-100 \
    --time 3:0:0 --job-name m3gnet-wbm-relax-IS2RE --mem 12000 \
    --output models/m3gnet/slurm_logs/slurm-%A-%a.out \
    --wrap "TF_CPP_MIN_LOG_LEVEL=2 python models/m3gnet/slurm_array_m3gnet_relax_wbm.py"
```

--time 2h is probably enough but missing indices are annoying so best be safe.

TF_CPP_MIN_LOG_LEVEL=2 means INFO and WARNING logs are not printed
https://stackoverflow.com/a/40982782

Requires M3GNet installation: pip install m3gnet
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

task_type = "IS2RE"
# task_type = "RS2RE"

job_id = os.environ.get("SLURM_JOB_ID", "debug")
job_array_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
# set large fallback job array size for fast testing/debugging
job_array_size = int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 10_000))

print(f"Job started running {datetime.now():%Y-%m-%d@%H-%M}")
print(f"{job_id = }")
print(f"{job_array_id = }")
print(f"{version('m3gnet') = }")

today = f"{datetime.now():%Y-%m-%d}"
out_dir = f"{ROOT}/data/{today}-m3gnet-wbm-relax-{task_type}"
os.makedirs(out_dir, exist_ok=True)
json_out_path = f"{out_dir}/{job_array_id}.json.gz"

if os.path.isfile(json_out_path):
    raise SystemExit(f"{json_out_path = } already exists, exciting early")

warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")
warnings.filterwarnings(action="ignore", category=UserWarning, module="tensorflow")


# %%
data_path = f"{ROOT}/data/2022-06-26-wbm-cses-and-initial-structures.json.gz"
print(f"Loading from {data_path=}")
df_wbm = pd.read_json(data_path).set_index("material_id")

df_this_job = np.array_split(df_wbm, job_array_size)[job_array_id]

run_params = dict(
    m3gnet_version=version("m3gnet"),
    job_id=job_id,
    job_array_id=job_array_id,
    data_path=data_path,
)
if wandb.run is None:
    wandb.login()

wandb.init(
    project="m3gnet",
    name=f"m3gnet-wbm-relax-{task_type}-{job_id}-{job_array_id}",
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


for material_id, structure in structures.items():
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
