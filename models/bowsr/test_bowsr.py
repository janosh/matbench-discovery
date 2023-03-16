# %%
from __future__ import annotations

import contextlib
import os
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
import wandb
from maml.apps.bowsr.model.megnet import MEGNet
from maml.apps.bowsr.optimizer import BayesianOptimizer
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import DEBUG, timestamp, today
from matbench_discovery.data import DATA_FILES, as_dict_handler
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

"""
To slurm submit this file: python path/to/file.py slurm-submit
Requires MEGNet and MAML installation: pip install megnet maml
https://github.com/materialsvirtuallab/maml
"""

task_type = "IS2RE"  # "RS2RE"
module_dir = os.path.dirname(__file__)
# set large job array size for smaller data splits and faster testing/debugging
slurm_array_task_count = 500
# see https://stackoverflow.com/a/55431306 for how to change array throttling
# post submission
slurm_max_parallel = 100
energy_model = "megnet"
job_name = f"bowsr-{energy_model}-wbm-{task_type}{'-debug' if DEBUG else ''}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")

data_path = {
    "IS2RE": DATA_FILES.wbm_initial_structures,
    "RS2RE": DATA_FILES.wbm_computed_structure_entries,
}[task_type]


slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    partition="skylake",
    account="LEE-SL3-CPU",
    time="12:0:0",
    # --time 2h is probably enough but best be safe.
    array=f"1-{slurm_array_task_count}%{slurm_max_parallel}",
    # --mem 12000 avoids slurmstepd: error: Detected 1 oom-kill event(s)
    # Some of your processes may have been killed by the cgroup out-of-memory handler.
    slurm_flags=("--mem", str(12_000)),
    # TF_CPP_MIN_LOG_LEVEL=2 means INFO and WARNING logs are not printed
    # https://stackoverflow.com/a/40982782
    pre_cmd="TF_CPP_MIN_LOG_LEVEL=2",
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", 0))
out_path = f"{out_dir}/bowsr-preds-{slurm_array_task_id}.json.gz"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")

print(f"\nJob started running {timestamp}")
print(f"{data_path = }")
print(f"{out_path = }")

df_in: pd.DataFrame = np.array_split(
    pd.read_json(data_path).set_index("material_id"), slurm_array_task_count
)[slurm_array_task_id - 1]


# %%
bayes_optim_kwargs = dict(
    relax_coords=True,
    relax_lattice=True,
    use_symmetry=True,
    seed=42,
)
optimize_kwargs = dict(n_init=100, n_iter=100, alpha=0.026**2)

run_params = dict(
    bayes_optim_kwargs=bayes_optim_kwargs,
    data_path=data_path,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    energy_model=energy_model,
    **{f"{dep}_version": version(dep) for dep in ("maml", "numpy", energy_model)},
    optimize_kwargs=optimize_kwargs,
    task_type=task_type,
    slurm_vars=slurm_vars,
)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %%
model = MEGNet()
relax_results: dict[str, dict[str, Any]] = {}
input_col = {"IS2RE": "initial_structure", "RS2RE": "relaxed_structure"}[task_type]

if task_type == "RS2RE":
    df_in[input_col] = [x["structure"] for x in df_in.computed_structure_entry]

structures = df_in[input_col].map(Structure.from_dict).to_dict()

for material_id in tqdm(structures, desc="Main loop", disable=None):
    structure = structures[material_id]
    if material_id in relax_results:
        continue
    try:
        optimizer = BayesianOptimizer(
            model=model, structure=structure, **bayes_optim_kwargs
        )
        optimizer.set_bounds()
        # reason for /dev/null: https://github.com/materialsvirtuallab/maml/issues/469
        with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull):
            optimizer.optimize(**optimize_kwargs)

        try:
            struct_bowsr, energy_bowsr = optimizer.get_optimized_structure_and_energy()
        except Exception as error:
            print(f"Failed to relax {material_id}: {error}")

        results = {
            f"e_form_per_atom_bowsr_{energy_model}": model.predict_energy(struct_bowsr),
            "structure_bowsr": struct_bowsr,
            f"energy_bowsr_{energy_model}": energy_bowsr,
        }

        relax_results[material_id] = results
    except Exception as exc:
        print(f"{material_id=} raised {exc=}")


# %%
df_out = pd.DataFrame(relax_results).T
df_out.index.name = "material_id"

df_out.reset_index().to_json(out_path, default_handler=as_dict_handler)

wandb.log_artifact(out_path, type=job_name)
