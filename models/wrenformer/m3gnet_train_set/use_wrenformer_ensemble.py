# %%
from __future__ import annotations

import os
from datetime import datetime

import pandas as pd
import wandb
from aviary.deploy import predict_from_wandb_checkpoints
from aviary.wrenformer.model import Wrenformer

from matbench_discovery.slurm import slurm_submit_python

__author__ = "Janosh Riebesell"
__date__ = "2022-09-05"

"""
Script that downloads checkpoints for an ensemble of Wrenformer models trained on the
per-atom energies in the M3GNet training set
https://figshare.com/articles/dataset/MPF_2021_2_8/19470599, then makes predictions on
WBM dataset, prints ensemble metrics and stores predictions to CSV.
"""

module_dir = os.path.dirname(__file__)
today = f"{datetime.now():%Y-%m-%d}"
ensemble_id = "wrenformer-m3gnet-trainset-ensemble-1"
run_name = f"{today}-{ensemble_id}-IS2RE"

slurm_submit_python(
    job_name=run_name,
    partition="ampere",
    account="LEE-SL3-GPU",
    time="1:0:0",
    log_dir=module_dir,
    slurm_flags=("--nodes", "1", "--gpus-per-node", "1"),
)


# %%
# download wbm-steps-summary.csv (23.31 MB)
data_path = "https://figshare.com/files/37570234?private_link=ff0ad14505f9624f0c05"
df = pd.read_csv(data_path).set_index("material_id")


target_col = "energy_per_atom"
df[target_col] = df.energy / df.n_sites

wandb.login()
wandb_api = wandb.Api()
runs = wandb_api.runs(
    "janosh/matbench-discovery", filters={"tags": {"$in": [ensemble_id]}}
)

assert len(runs) == 10, f"Expected 10 runs, got {len(runs)} for {ensemble_id=}"

df, ensemble_metrics = predict_from_wandb_checkpoints(
    runs, df, input_col="wyckoff", target_col=target_col, model_cls=Wrenformer
)

df.round(6).to_csv(f"{module_dir}/{today}-{run_name}-preds.csv")
