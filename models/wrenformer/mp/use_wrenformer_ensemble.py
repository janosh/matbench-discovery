# %%
from __future__ import annotations

import os
from datetime import datetime

import pandas as pd
import wandb
from aviary.deploy import predict_from_wandb_checkpoints
from aviary.wrenformer.data import df_to_in_mem_dataloader
from aviary.wrenformer.model import Wrenformer

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

"""
Script that downloads checkpoints for an ensemble of Wrenformer models trained on the MP
formation energies, then makes predictions on some dataset, prints ensemble metrics and
stores predictions to CSV.
"""

module_dir = os.path.dirname(__file__)
today = f"{datetime.now():%Y-%m-%d}"


# %%
# download wbm-steps-summary.csv (23.31 MB)
data_path = "https://figshare.com/files/37570234?private_link=ff0ad14505f9624f0c05"
df = pd.read_csv(data_path).set_index("material_id")

target_col = "e_form_per_atom"
input_col = "wyckoff"
assert target_col in df, f"{target_col=} not in {list(df)}"
assert input_col in df, f"{input_col=} not in {list(df)}"

wandb.login()
wandb_api = wandb.Api()
ensemble_id = "wrenformer-e_form-ensemble-1"
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
    runs, data_loader, df=df, model_class=Wrenformer
)

df.round(6).to_csv(f"{module_dir}/{today}-{ensemble_id}-preds-{target_col}.csv")
