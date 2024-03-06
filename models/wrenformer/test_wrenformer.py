"""Download WandB checkpoints for an ensemble of Wrenformer models trained on all MP
formation energies, then makes predictions on some dataset, prints ensemble metrics and
saves predictions to CSV.
"""

# %%
from __future__ import annotations

import os
import sys
from importlib.metadata import version

import torch
import wandb
from aviary.predict import predict_from_wandb_checkpoints
from aviary.wrenformer.data import df_to_in_mem_dataloader
from aviary.wrenformer.model import Wrenformer

from matbench_discovery import CHECKPOINT_DIR, WANDB_PATH, today
from matbench_discovery.data import df_wbm
from matbench_discovery.enums import Key, Task
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"


task_type = Task.IS2RE
debug = "slurm-submit" in sys.argv
job_name = f"test-wrenformer-wbm-{task_type}"
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    partition="ampere",
    account="LEE-SL3-GPU",
    time="2:0:0",
    out_dir=out_dir,
    slurm_flags="--nodes 1 --gpus-per-node 1",
)


# %%
df = df_wbm.dropna(subset=Key.init_wyckoff)

assert Key.e_form in df, f"{Key.e_form=} not in {list(df)}"
assert Key.wyckoff in df, f"{Key.wyckoff=} not in {list(df)}"


# %%
filters = {
    "created_at": {"$gt": "2022-11-15", "$lt": "2022-11-16"},
    "display_name": {"$regex": "wrenformer-"},
}
runs = wandb.Api().runs(WANDB_PATH, filters=filters)
expected_runs = 10
assert (
    len(runs) == expected_runs
), f"{expected_runs=}, got {len(runs)} filtering {WANDB_PATH=} with {filters=}"

for idx, run in enumerate(runs):
    for key, val in run.config.items():
        if val == runs[0].config[key] or key.startswith(("slurm_", "timestamp")):
            continue
        raise ValueError(
            f"Run configs not identical: runs[{idx}][{key}]={val}, {runs[0][key]=}"
        )

run_params = dict(
    df=dict(shape=str(df.shape), columns=", ".join(df)),
    versions={dep: version(dep) for dep in ("aviary", "numpy", "torch")},
    ensemble_size=len(runs),
    task_type=task_type,
    target_col=Key.e_form,
    input_col=Key.wyckoff,
    wandb_run_filters=filters,
    slurm_vars=slurm_vars,
    training_run_ids=[run.id for run in runs],
)

try:  # load checkpoint to get number of parameters
    runs[0].file("checkpoint.pth").download(root=module_dir)
    state_dict = torch.load(f"{module_dir}/checkpoint.pth", map_location="cpu")
    model = Wrenformer(**state_dict["model_params"])
    run_params[Key.model_params] = model.num_params
except Exception as exc:
    print(exc)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %%
data_loader = df_to_in_mem_dataloader(
    df=df,
    cache_dir=CHECKPOINT_DIR,
    target_col=Key.e_form,
    batch_size=1024,
    input_col=Key.wyckoff,
    embedding_type="wyckoff",
    shuffle=False,  # False is default but best be explicit
)

df, ensemble_metrics = predict_from_wandb_checkpoints(
    runs,
    data_loader=data_loader,
    df=df,
    model_cls=Wrenformer,
    target_col=Key.e_form,
)
df = df.round(4)

slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")
df.to_csv(f"{out_dir}/{job_name}-preds-{slurm_array_job_id}.csv.gz")


# %%
pred_col = f"{Key.e_form}_pred_ens"
assert pred_col in df, f"{pred_col=} not in {list(df)}"
table = wandb.Table(dataframe=df[[Key.e_form, pred_col]].reset_index())


# %%
MAE = ensemble_metrics.MAE.mean()
R2 = ensemble_metrics.R2.mean()

title = f"Wrenformer {task_type} ensemble={len(runs)} {MAE=:.4} {R2=:.4}"

wandb_scatter(table, fields=dict(x=Key.e_form, y=pred_col), title=title)
