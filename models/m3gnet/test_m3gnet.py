"""Get M3GNet formation energy predictions on WBM test set.
To slurm submit this file: python path/to/file.py slurm-submit
Requires M3GNet installation: pip install m3gnet
https://github.com/materialsvirtuallab/m3gnet.
"""


# %%
from __future__ import annotations

import os
import warnings
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
import wandb
from m3gnet.models import Relaxer
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import DEBUG, timestamp, today
from matbench_discovery.data import DATA_FILES, as_dict_handler
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

task_type = "IS2RE"  # "RS2RE"
module_dir = os.path.dirname(__file__)
# set large job array size for fast testing/debugging
slurm_array_task_count = 100
job_name = f"m3gnet-wbm-{task_type}{'-debug' if DEBUG else ''}"
out_dir = os.environ.get("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time="3:0:0",
    array=f"1-{slurm_array_task_count}",
    slurm_flags=("--mem", str(12_000)),
    # TF_CPP_MIN_LOG_LEVEL=2 means INFO and WARNING logs are not printed
    # https://stackoverflow.com/a/40982782
    pre_cmd="TF_CPP_MIN_LOG_LEVEL=2",
)


# %%
slurm_array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 3))
out_path = f"{out_dir}/m3gnet-preds-{slurm_array_task_id}.json.gz"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")

warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")
warnings.filterwarnings(action="ignore", category=UserWarning, module="tensorflow")


# %%
data_path = {
    "IS2RE": DATA_FILES.wbm_initial_structures,
    "RS2RE": DATA_FILES.wbm_computed_structure_entries,
}[task_type]
print(f"\nJob started running {timestamp}")
print(f"{data_path=}")
df_wbm = pd.read_json(data_path).set_index("material_id")

df_in: pd.DataFrame = np.array_split(df_wbm, slurm_array_task_count)[
    slurm_array_task_id - 1
]

run_params = dict(
    data_path=data_path,
    m3gnet_version=version("m3gnet"),
    numpy_version=version("numpy"),
    task_type=task_type,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    slurm_vars=slurm_vars,
)

run_name = f"{job_name}-{slurm_array_task_id}"
wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %%
megnet = Relaxer()  # load default pre-trained M3GNet model
relax_results: dict[str, dict[str, Any]] = {}
input_col = {"IS2RE": "initial_structure", "RS2RE": "relaxed_structure"}[task_type]

if task_type == "RS2RE":
    df_in[input_col] = [x["structure"] for x in df_in.computed_structure_entry]

structures = df_in[input_col].map(Structure.from_dict).to_dict()

for material_id in tqdm(structures, disable=None):
    if material_id in relax_results:
        continue
    try:
        relax_result = megnet.relax(structures[material_id])
    except Exception as error:
        print(f"Failed to relax {material_id}: {error}")
        continue
    relax_dict = {
        "m3gnet_structure": relax_result["final_structure"],
        "m3gnet_trajectory": relax_result["trajectory"].__dict__,
    }

    relax_results[material_id] = relax_dict


# %%
df_out = pd.DataFrame(relax_results).T
df_out.index.name = "material_id"

df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)

wandb.log_artifact(out_path, type=f"m3gnet-wbm-{task_type}")
