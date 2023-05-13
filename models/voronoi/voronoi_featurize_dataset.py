"""Featurize MP training and WBM test structures with Magpie composition-based and
Voronoi tessellation structure-based features.
"""


# %%
import os
import sys
import warnings
from importlib.metadata import version

import numpy as np
import pandas as pd
import wandb
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import DEBUG, today
from matbench_discovery.data import DATA_FILES
from matbench_discovery.slurm import slurm_submit
from models.voronoi import featurizer

__author__ = "Janosh Riebesell"
__date__ = "2022-10-31"


data_name = "mp"
data_path = {
    "wbm": DATA_FILES.wbm_initial_structures,
    "mp": DATA_FILES.mp_computed_structure_entries,
}[data_name]

input_col = "initial_structure"
# input_col = "relaxed_structure"
debug = "slurm-submit" in sys.argv
job_name = f"voronoi-features-{data_name}{'-debug' if DEBUG else ''}"
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")
slurm_array_task_count = 50


slurm_vars = slurm_submit(
    job_name=job_name,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time="12:0:0",
    array=f"1-{slurm_array_task_count}",
    slurm_flags=("--mem", "15G") if data_name == "mp" else (),
    out_dir=out_dir,
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
run_name = f"{job_name}-{slurm_array_task_id}"
out_path = f"{out_dir}/{run_name}.csv.bz2"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")

print(f"{data_path=}")
df_in: pd.DataFrame = np.array_split(
    pd.read_json(data_path).set_index("material_id"), slurm_array_task_count
)[slurm_array_task_id - 1]

if data_name == "mp":  # extract structure dicts from ComputedStructureEntry
    struct_dicts = [x["structure"] for x in df_in.entry]
elif data_name == "wbm" and input_col == "relaxed_structure":
    struct_dicts = [x["structure"] for x in df_in.computed_structure_entry]
elif data_name == "wbm" and input_col == "initial_structure":
    struct_dicts = df_in.initial_structure

df_in[input_col] = [Structure.from_dict(x) for x in tqdm(struct_dicts, disable=None)]


# %%
run_params = dict(
    data_path=data_path,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    input_col=input_col,
    slurm_vars=slurm_vars,
    out_path=out_path,
    **{f"{dep}_version": version(dep) for dep in ("matminer", "numpy")},
)

wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %% prints lots of pymatgen warnings
# > No electronegativity for Ne. Setting to NaN. This has no physical meaning, ...
warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")

df_features = featurizer.featurize_dataframe(df_in, input_col, ignore_errors=True)[
    featurizer.feature_labels()
].round(4)


# %%
df_features.to_csv(out_path)

wandb.log({"voronoi_features": wandb.Table(dataframe=df_features)})
