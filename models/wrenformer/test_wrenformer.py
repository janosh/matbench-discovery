# %%
from __future__ import annotations

import os
import sys
from importlib.metadata import version

import pandas as pd
import wandb
from aviary.predict import predict_from_wandb_checkpoints
from aviary.wrenformer.data import df_to_in_mem_dataloader
from aviary.wrenformer.model import Wrenformer

from matbench_discovery import CHECKPOINT_DIR, DEBUG, ROOT, WANDB_PATH, today
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

"""
Download WandB checkpoints for an ensemble of Wrenformer models trained on all MP
formation energies, then makes predictions on some dataset, prints ensemble metrics and
saves predictions to CSV.
"""

task_type = "IS2RE"
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-summary.csv"
debug = "slurm-submit" in sys.argv
job_name = f"test-wrenformer-wbm-{task_type}{'-debug' if DEBUG else ''}"
module_dir = os.path.dirname(__file__)
out_dir = os.environ.get("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    partition="ampere",
    account="LEE-SL3-GPU",
    time="2:0:0",
    out_dir=out_dir,
    slurm_flags="--nodes 1 --gpus-per-node 1",
)


# %%
target_col = "e_form_per_atom_mp2020_corrected"
input_col = "wyckoff_spglib"
df = pd.read_csv(data_path).dropna(subset=input_col).set_index("material_id")

assert target_col in df, f"{target_col=} not in {list(df)}"
assert input_col in df, f"{input_col=} not in {list(df)}"


# %%
filters = {
    "created_at": {"$gt": "2022-11-15", "$lt": "2022-11-16"},
    "display_name": {"$regex": "wrenformer-robust"},
}
runs = wandb.Api().runs(WANDB_PATH, filters=filters)

assert len(runs) == 10, f"Expected 10 runs, got {len(runs)} for {filters=}"
for idx, run in enumerate(runs):
    for key, val in run.config.items():
        if val == runs[0].config[key] or key.startswith(("slurm_", "timestamp")):
            continue
        raise ValueError(
            f"Run configs not identical: runs[{idx}][{key}]={val}, {runs[0][key]=}"
        )

run_params = dict(
    data_path=data_path,
    df=dict(shape=str(df.shape), columns=", ".join(df)),
    aviary_version=version("aviary"),
    numpy_version=version("numpy"),
    torch_version=version("torch"),
    ensemble_size=len(runs),
    task_type=task_type,
    target_col=target_col,
    input_col=input_col,
    wandb_run_filters=filters,
    slurm_vars=slurm_vars,
)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %%
data_loader = df_to_in_mem_dataloader(
    df=df,
    cache_dir=CHECKPOINT_DIR,
    target_col=target_col,
    batch_size=1024,
    input_col=input_col,
    embedding_type="wyckoff",
    shuffle=False,  # False is default but best be explicit
)

df, ensemble_metrics = predict_from_wandb_checkpoints(
    runs, data_loader=data_loader, df=df, model_cls=Wrenformer, target_col=target_col
)
df = df.round(4)

slurm_job_id = os.environ.get("SLURM_JOB_ID", "debug")
df.to_csv(f"{out_dir}/{job_name}-preds-{slurm_job_id}.csv")


# %%
pred_col = f"{target_col}_pred_ens"
assert pred_col in df, f"{pred_col=} not in {list(df)}"
table = wandb.Table(dataframe=df[[target_col, pred_col]].reset_index())


# %%
MAE = ensemble_metrics.MAE.mean()
R2 = ensemble_metrics.R2.mean()

title = f"Wrenformer {task_type} ensemble={len(runs)} {MAE=:.4} {R2=:.4}"

wandb_scatter(table, fields=dict(x=target_col, y=pred_col), title=title)
