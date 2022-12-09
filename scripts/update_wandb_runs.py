# %%
import pandas as pd
import wandb
from wandb.wandb_run import Run

from matbench_discovery import WANDB_PATH

"""
Update run metadata recorded on Weights and Biases
https://wandb.ai/janosh/matbench-discovery.
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-09-21"


# %%
filters = dict(display_name={"$regex": "voronoi-featurize"})
runs = wandb.Api().runs(WANDB_PATH, filters=filters)

print(f"matching runs: {len(runs)}")


# %%
df = pd.DataFrame([run.config | dict(run.summary) for run in runs])
df["display_name"] = [run.display_name for run in runs]


# %%
df.isna().sum()


# %% --- Update run metadata ---
updated_runs: list[Run] = []
wet_run = input("Wet run or dry run? [w/d] ").lower().startswith("w")

for idx, run in enumerate(runs, 1):
    old_config, new_config = run.config.copy(), run.config.copy()

    new_display_name = run.display_name.replace("featurize", "features")

    for x in ("IS2RE", "ES2RE"):
        if x in run.display_name:
            new_config["task_type"] = x

    if "SLURM_JOB_ID" in new_config:
        new_config["slurm_job_id"] = new_config.pop("SLURM_JOB_ID")

    if "SLURM_ARRAY_TASK_ID" in new_config:
        new_config["slurm_array_task_id"] = new_config.pop("SLURM_ARRAY_TASK_ID")

    if old_config != new_config or new_display_name != run.display_name:
        print(f"\nrun {idx}/{len(runs)}: {run.display_name}")

        if new_display_name != run.display_name:
            print(f"{new_display_name=}")

        for key in set(old_config) | set(new_config):
            old_val = old_config.get(key)
            new_val = new_config.get(key)
            if new_val != old_val:
                print(f"{key}: {old_val} => {new_val}")

        updated_runs.append(run)

        if wet_run:
            run.display_name = new_display_name
            run.config = new_config
            run.update()

print(f"\n{'' if wet_run else 'dry run: would have'} updated {len(updated_runs)} runs")
