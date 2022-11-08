# %%
from __future__ import annotations

import os
from datetime import datetime

import pandas as pd
import wandb
from aviary.deploy import predict_from_wandb_checkpoints
from aviary.wrenformer.data import df_to_in_mem_dataloader
from aviary.wrenformer.model import Wrenformer

from matbench_discovery import ROOT
from matbench_discovery.slurm import slurm_submit_python

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

"""
Script that downloads checkpoints for an ensemble of Wrenformer models trained on the MP
formation energies, then makes predictions on some dataset, prints ensemble metrics and
stores predictions to CSV.
"""

module_dir = os.path.dirname(__file__)
today = f"{datetime.now():%Y-%m-%d}"
ensemble_id = "wrenformer-e_form-ensemble-1"
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
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-summary.csv"
target_col = "e_form_per_atom_mp2020_corrected"
input_col = "wyckoff_spglib"
df = pd.read_csv(data_path).dropna(subset=input_col).set_index("material_id")

assert target_col in df, f"{target_col=} not in {list(df)}"
assert input_col in df, f"{input_col=} not in {list(df)}"

wandb.login()
wandb_api = wandb.Api()
runs = wandb_api.runs(
    "janosh/matbench-discovery", filters={"tags": {"$in": [ensemble_id]}}
)

assert len(runs) == 10, f"Expected 10 runs, got {len(runs)} for {ensemble_id=}"

data_loader = df_to_in_mem_dataloader(
    df=df,
    target_col=target_col,
    batch_size=1024,
    input_col=input_col,
    embedding_type="wyckoff",
    shuffle=False,  # False is default but best be explicit
)

df, ensemble_metrics = predict_from_wandb_checkpoints(
    runs, data_loader=data_loader, df=df, model_cls=Wrenformer
)

df.round(6).to_csv(f"{module_dir}/{today}-{run_name}-preds.csv")
