# %%
from __future__ import annotations

import os
from datetime import datetime

import pandas as pd
import wandb
from aviary.wrenformer.deploy import deploy_wandb_checkpoints

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


# %%
# download wbm-steps-summary.csv (23.31 MB)
data_path = "https://figshare.com/files/37570234?private_link=ff0ad14505f9624f0c05"
df = pd.read_csv(data_path).set_index("material_id")


target_col = "energy_per_atom"
df[target_col] = df.energy / df.n_sites

wandb.login()
wandb_api = wandb.Api()
ensemble_id = "wrenformer-m3gnet-trainset-ensemble-1"
runs = wandb_api.runs(
    "aviary/wrenformer-on-m3gnet-trainset", filters={"tags": {"$in": [ensemble_id]}}
)

assert len(runs) == 10, f"Expected 10 runs, got {len(runs)} for {ensemble_id=}"

df, ensemble_metrics = deploy_wandb_checkpoints(
    runs, df, input_col="wyckoff", target_col=target_col
)

df.round(6).to_csv(f"{module_dir}/{today}-{ensemble_id}-preds-{target_col}.csv")
