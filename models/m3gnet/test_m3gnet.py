"""Get M3GNet formation energy predictions on WBM test set.
To slurm submit this file: python path/to/file.py slurm-submit
Requires M3GNet installation: pip install m3gnet
https://github.com/materialsvirtuallab/m3gnet.
"""

# %%
import os
import warnings
from importlib.metadata import version
from typing import Any, Literal

import numpy as np
import pandas as pd
import wandb
from m3gnet.models import Relaxer
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import ROOT, timestamp, today
from matbench_discovery.data import DataFiles, as_dict_handler
from matbench_discovery.enums import Task
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

task_type = Task.IS2RE
module_dir = os.path.dirname(__file__)
model_name = "m3gnet"
# direct: DIRECT cluster sampling, ms: manual sampling
model_type: Literal["orig", "direct", "manual-sampling"] = "orig"
# set large job array size for smaller data splits and faster testing/debugging
slurm_array_task_count = 50
record_traj = False
job_name = f"{model_name}-{model_type}-wbm-{task_type}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    account="matgen",
    time="11:55:0",
    array=f"1-{slurm_array_task_count}",
    slurm_flags="--qos shared --constraint cpu --mem 16G",
    # TF_CPP_MIN_LOG_LEVEL=2 means INFO and WARNING logs are not printed
    # https://stackoverflow.com/a/40982782
    pre_cmd="TF_CPP_MIN_LOG_LEVEL=2",
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "3"))
slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")
out_path = f"{out_dir}/{slurm_array_job_id}-{slurm_array_task_id:>03}.json.gz"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exciting early")

warnings.filterwarnings(action="ignore", category=UserWarning, module="tensorflow")


# %%
data_path = {
    Task.IS2RE: DataFiles.wbm_initial_structures.path,
    Task.RS2RE: DataFiles.wbm_computed_structure_entries.path,
}[task_type]
print(f"\nJob {job_name} started {timestamp}")
print(f"{data_path=}")
e_pred_col = f"{model_name}_{model_type}_energy"

df_in = pd.read_json(data_path).set_index(Key.mat_id)
if slurm_array_task_count > 1:
    df_in = np.array_split(df_in, slurm_array_task_count)[slurm_array_task_id - 1]

checkpoint = None
if model_type == "direct":
    checkpoint = f"{ROOT}/models/{model_name}/2023-05-26-DI-DFTstrictF10-TTRS-128U-442E"
if model_type == "ms":
    checkpoint = f"{ROOT}/models/{model_name}/2023-05-26-MS-DFTstrictF10-128U-154E"
relax_results: dict[str, dict[str, Any]] = {}
m3gnet = Relaxer(potential=checkpoint)  # load pre-trained M3GNet model

run_params = {
    "data_path": data_path,
    "versions": {dep: version(dep) for dep in ("m3gnet", "numpy")},
    Key.task_type: task_type,
    "df": {"shape": str(df_in.shape), "columns": ", ".join(df_in)},
    "slurm_vars": slurm_vars,
    Key.model_params: sum(
        np.prod(weight.shape) for weight in m3gnet.potential.model.trainable_weights
    ),
    "checkpoint": checkpoint,
    "model_type": model_type,
    "out_path": out_path,
    "job_name": job_name,
    "record_traj": record_traj,
}

run_name = f"{job_name}-{slurm_array_task_id}"
wandb.init(project="matbench-discovery", name=run_name, config=run_params)


# %%
input_col = {Task.IS2RE: Key.init_struct, Task.RS2RE: Key.final_struct}[task_type]

if task_type == Task.RS2RE:
    df_in[input_col] = [cse["structure"] for cse in df_in[Key.cse]]

structures = df_in[input_col].map(Structure.from_dict).to_dict()

for material_id in tqdm(structures, desc="Relaxing"):
    if material_id in relax_results:
        continue
    try:
        result = m3gnet.relax(structures[material_id])
        relax_results[material_id] = {
            f"{model_name}_{model_type}_structure": result["final_structure"],
            e_pred_col: result["trajectory"].energies[-1],
        }
        if record_traj:
            relax_results[f"{model_name}_{model_type}_trajectory"] = result[
                "trajectory"
            ].__dict__
    except Exception as exc:
        print(f"Failed to relax {material_id}: {exc!r}")


# %%
df_out = pd.DataFrame(relax_results).T
df_out.index.name = Key.mat_id

df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)

wandb.log_artifact(out_path, type=f"{model_name}-wbm-{task_type}")
