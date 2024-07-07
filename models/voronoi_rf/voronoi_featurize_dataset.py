"""Featurize MP training and WBM test structures with Magpie composition-based and
Voronoi tessellation structure-based features.
"""

# %%
import os
import sys
from importlib.metadata import version

import numpy as np
import pandas as pd
import wandb
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import ROOT, today
from matbench_discovery.data import DATA_FILES
from matbench_discovery.slurm import slurm_submit

sys.path.append(f"{ROOT}/models")

from voronoi_rf import featurizer

__author__ = "Janosh Riebesell"
__date__ = "2022-10-31"


# %%
data_name = "mp"
data_path = {
    "wbm": DATA_FILES.wbm_initial_structures,
    "mp": DATA_FILES.mp_computed_structure_entries,
}[data_name]

input_col = Key.init_struct  # or Key.final_struct
debug = "slurm-submit" in sys.argv
job_name = f"voronoi-features-{data_name}"
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")
slurm_array_task_count = 50


slurm_vars = slurm_submit(
    job_name=job_name,
    account="matgen",
    time="11:55:0",
    array=f"1-{slurm_array_task_count}",
    slurm_flags=("--mem", "15G") if data_name == "mp" else (),
    out_dir=out_dir,
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
run_name = f"{job_name}-{slurm_array_task_id}"
out_path = f"{out_dir}/{run_name}.csv.bz2"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exciting early")

print(f"{data_path=}")
df_in = pd.read_json(data_path).set_index(Key.mat_id)
if slurm_array_task_count > 1:
    df_in = np.array_split(df_in, slurm_array_task_count)[slurm_array_task_id - 1]

if data_name == "mp":  # extract structure dicts from ComputedStructureEntry
    struct_dicts = [cse["structure"] for cse in df_in.entry]
elif data_name == "wbm" and input_col == Key.final_struct:
    struct_dicts = [cse["structure"] for cse in df_in[Key.cse]]
elif data_name == "wbm" and input_col == Key.init_struct:
    struct_dicts = df_in[Key.init_struct]
else:
    raise ValueError(f"Invalid {data_name=}, {input_col=} combo")

df_in[input_col] = [
    Structure.from_dict(dct) for dct in tqdm(struct_dicts, disable=None)
]


# %%
run_params = dict(
    data_path=data_path,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    input_col=input_col,
    slurm_vars=slurm_vars,
    out_path=out_path,
    versions={dep: version(dep) for dep in ("matminer", "numpy")},
)

wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %%
df_features = featurizer.featurize_dataframe(df_in, input_col, ignore_errors=True)[
    featurizer.feature_labels()
].round(4)


# %%
df_features.to_csv(out_path)

wandb.log({"voronoi_features": wandb.Table(dataframe=df_features)})
