# %%
from __future__ import annotations

import os
from datetime import datetime
from importlib.metadata import version

import pandas as pd
import wandb
from aviary.deploy import predict_from_wandb_checkpoints
from aviary.wrenformer.data import df_to_in_mem_dataloader
from aviary.wrenformer.model import Wrenformer

from matbench_discovery import ROOT
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

"""
Download wandb checkpoints for an ensemble of Wrenformer models trained on MP
formation energies, then makes predictions on some dataset, prints ensemble metrics and
stores predictions to CSV.
"""

module_dir = os.path.dirname(__file__)
today = f"{datetime.now():%Y-%m-%d}"
task_type = "IS2RE"
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-summary.csv"
job_name = "wrenformer-wbm-IS2RE"

slurm_vars = slurm_submit(
    job_name=job_name,
    partition="ampere",
    account="LEE-SL3-GPU",
    time="2:0:0",
    log_dir=module_dir,
    slurm_flags=("--nodes", "1", "--gpus-per-node", "1"),
)


# %%
target_col = "e_form_per_atom_mp2020_corrected"
input_col = "wyckoff_spglib"
df = pd.read_csv(data_path).dropna(subset=input_col).set_index("material_id")

assert target_col in df, f"{target_col=} not in {list(df)}"
assert input_col in df, f"{input_col=} not in {list(df)}"


# %%
wandb.login()

filters = {
    "created_at": {"$gt": "2022-11-15", "$lt": "2022-11-16"},
    "display_name": {"$regex": "wrenformer-robust"},
}
runs = wandb.Api().runs("janosh/matbench-discovery", filters=filters)

assert len(runs) == 10, f"Expected 10 runs, got {len(runs)} for {filters=}"
for idx, run in enumerate(runs):
    for key, val in run.config.items():
        if val == runs[0][key] or key.startswith(("slurm_", "timestamp")):
            continue
        raise ValueError(
            f"Configs not identical: runs[{idx}][{key}]={val}, {runs[0][key]=}"
        )

run_params = dict(
    data_path=data_path,
    df=dict(shape=str(df.shape), columns=", ".join(df)),
    aviary_version=version("aviary"),
    ensemble_size=len(runs),
    task_type=task_type,
    target_col=target_col,
    input_col=input_col,
    filters=filters,
    slurm_vars=slurm_vars,
)

slurm_job_id = os.environ.get("SLURM_JOB_ID", "debug")
wandb.init(
    project="matbench-discovery", name=f"{job_name}-{slurm_job_id}", config=run_params
)


# %%
data_loader = df_to_in_mem_dataloader(
    df=df,
    target_col=target_col,
    batch_size=1024,
    input_col=input_col,
    embedding_type="wyckoff",
    shuffle=False,  # False is default but best be explicit
)

df, ensemble_metrics = predict_from_wandb_checkpoints(
    runs, data_loader=data_loader, df=df, model_cls=Wrenformer, target_col=target_col
)

df.to_csv(f"{module_dir}/{today}-{job_name}-preds.csv")


# %%
pred_col = f"{target_col}_pred_ens"
assert pred_col in df, f"{pred_col=} not in {list(df)}"
table = wandb.Table(dataframe=df[[target_col, pred_col]])


# %%
MAE = ensemble_metrics["MAE"]
R2 = ensemble_metrics["R2"]

title = rf"Wrenformer {task_type} ensemble={len(runs)} {MAE=:.4} {R2=:.4}"
fields = dict(x=target_col, y=pred_col, title=title)
print(title)

scatter_plot = wandb.plot_table(
    vega_spec_name="janosh/scatter-parity", data_table=table, fields=fields
)

wandb.log({"true_pred_scatter": scatter_plot})
