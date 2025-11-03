"""Test MatRIS model for matbench-discovery"""
# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = [
# "torch==2.6.0",
# "numpy==2.3.4",
# "ase==3.26.0",
# "pymatgen==2025.10.7",
# "pymatviz==0.17.2",
# "scikit-learn==1.7.2",
# "tqdm==4.67.1",
# "matbench-discovery"==1.3.1,
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../../", editable = true }
# ///
# %%
import os
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
import torch
from matris.applications.relax import StructOptimizer
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import timestamp, today
from matbench_discovery.data import as_dict_handler
from matbench_discovery.enums import DataFiles, Task
from matbench_discovery.hpc import slurm_submit

task_type = Task.IS2RE
module_dir = os.path.dirname(__file__)
slurm_array_task_count = 45

model_name = "MatRIS_10M_MP"  # or MatRIS_10M_OAM
device = "cuda" if torch.cuda.is_available() else "cpu"
matris = StructOptimizer(model=model_name, device=device)

job_name = f"{model_name}_Discovery/{today}-wbm-{task_type}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    time="47:55:0",
    partition="gpu",
    array=f"1-{slurm_array_task_count}",
    slurm_flags="--nodes 1 --gres gpu:1 --cpus-per-task 1",
    submit_as_temp_file=False,
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "1"))
slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")

out_path = f"{out_dir}/{slurm_array_job_id}-{slurm_array_task_id:>03}.json.gz"
if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exiting early")


# %%
data_path = {
    Task.RS2RE: DataFiles.wbm_computed_structure_entries.path,
    Task.IS2RE: DataFiles.wbm_initial_structures.path,
}[task_type]
print(f"\nJob {job_name} started {timestamp}")
print(f"{data_path=}")
e_pred_col = "matris_energy"
max_steps = 500
fmax = 0.02

df_in = pd.read_json(data_path, lines=True).set_index(Key.mat_id)
if slurm_array_task_count > 1:
    df_in = np.array_split(df_in, slurm_array_task_count)[slurm_array_task_id - 1]

run_params = {
    "data_path": data_path,
    "versions": {dep: version(dep) for dep in ("numpy", "torch")},
    Key.task_type: task_type,
    "df": {"shape": str(df_in.shape), "columns": ", ".join(df_in)},
    "slurm_vars": slurm_vars,
    "max_steps": max_steps,
    "fmax": fmax,
    "device": device,
    # Key.model_params: matris.n_params,
}

run_name = f"{job_name}-{slurm_array_task_id}"


# %%
relax_results: dict[str, dict[str, Any]] = {}
input_col = {Task.IS2RE: Key.initial_struct, Task.RS2RE: Key.final_struct}[task_type]
structures = df_in[input_col].map(Structure.from_dict).to_dict()

for material_id in tqdm(structures, desc="Relaxing"):
    if material_id in relax_results:
        continue
    try:
        relax_result = matris.relax(
            structures[material_id],
            verbose=False,
            steps=max_steps,
            fmax=fmax,
            relax_cell=max_steps > 0,
            ase_filter="FrechetCellFilter",
        )
        relax_results[material_id] = {
            e_pred_col: relax_result["trajectory"].energies[-1]
        }
        if max_steps > 0:
            relax_struct = relax_result["final_structure"]
            relax_results[material_id]["matris_structure"] = relax_struct
    except Exception as exc:
        filename = f"error_crystals/{model_name}/"
        if not os.path.exists(filename):
            os.makedirs(filename)
        structures[material_id].to(filename=f"{filename}{material_id}.cif")

        print(f"Failed to relax {material_id}: {exc!r}")


# %%
df_out = pd.DataFrame(relax_results).T
df_out.index.name = Key.mat_id

if max_steps == 0:
    df_out.add_suffix("_no_relax").to_csv(out_path.replace(".json.gz", ".csv.gz"))
else:
    df_out.reset_index().to_json(
        out_path, default_handler=as_dict_handler, orient="records", lines=True
    )
