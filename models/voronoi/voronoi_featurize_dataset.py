# %%
import os
import warnings
from datetime import datetime

import numpy as np
import pandas as pd
import wandb
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import ROOT, as_dict_handler
from matbench_discovery.slurm import slurm_submit
from models.voronoi import featurizer

today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)

data_name = "mp"  # "mp"
if data_name == "wbm":
    data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
    input_col = "initial_structure"
elif data_name == "mp":
    data_path = f"{ROOT}/data/mp/2022-09-16-mp-computed-structure-entries.json.gz"
    input_col = "structure"

slurm_array_task_count = 30
job_name = f"voronoi-features-{data_name}"
log_dir = f"{module_dir}/{today}-{job_name}"

slurm_vars = slurm_submit(
    job_name=job_name,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time=(slurm_max_job_time := "12:0:0"),
    array=f"1-{slurm_array_task_count}",
    slurm_flags=("--mem", "20G") if data_name == "mp" else (),
    log_dir=log_dir,
)


# %%
slurm_array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
out_path = f"{log_dir}/{job_name}.csv.bz2"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")

print(f"{data_path=}")
df = pd.read_json(data_path).set_index("material_id")
df_this_job: pd.DataFrame = np.array_split(df, slurm_array_task_count)[
    slurm_array_task_id - 1
]

if data_name == "mp":  # extract structure dicts from ComputedStructureEntry
    struct_dicts = [x["structure"] for x in df_this_job.entry]
if data_name == "wbm":
    struct_dicts = df_this_job.initial_structure

df_this_job[input_col] = [
    Structure.from_dict(x) for x in tqdm(struct_dicts, disable=None)
]


# %%
run_params = dict(
    data_path=data_path,
    df=dict(shape=str(df_this_job.shape), columns=", ".join(df_this_job)),
    input_col=input_col,
    slurm_vars=slurm_vars | dict(slurm_max_job_time=slurm_max_job_time),
)
if wandb.run is None:
    wandb.login()

slurm_job_id = os.environ.get("SLURM_JOB_ID", "debug")
wandb.init(
    project="matbench-discovery",
    name=f"{job_name}-{slurm_job_id}-{slurm_array_task_id}",
    config=run_params,
)


# %% prints lots of pymatgen warnings
# > No electronegativity for Ne. Setting to NaN. This has no physical meaning, ...
warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")

df_features = featurizer.featurize_dataframe(
    df_this_job, input_col, ignore_errors=True, pbar=dict(position=0, leave=True)
)


# %%
df_features[featurizer.feature_labels()].to_csv(
    out_path, default_handler=as_dict_handler
)

wandb.log({"voronoi_features": wandb.Table(dataframe=df_features)})
