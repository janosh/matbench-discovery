"""Download WandB checkpoints for an ensemble of Wrenformer models trained on all MP
formation energies, then makes predictions on some dataset, prints ensemble metrics and
saves predictions to CSV.
"""

# %%
import os
import sys
from importlib.metadata import version

import torch
import wandb
from aviary.predict import predict_from_wandb_checkpoints
from aviary.wrenformer.data import df_to_in_mem_dataloader
from aviary.wrenformer.model import Wrenformer
from pymatviz.enums import Key

from matbench_discovery import WANDB_PATH, today
from matbench_discovery.data import df_wbm
from matbench_discovery.enums import MbdKey, Model, Task
from matbench_discovery.hpc import slurm_submit
from matbench_discovery.plots import wandb_scatter

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"


task_type = Task.IS2RE
debug = "slurm-submit" in sys.argv
model_name = f"{Model.wrenformer}-ens=10"
job_name = f"{model_name}/{today}-wbm-{task_type}"
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    account="matgen",
    time="2:0:0",
    out_dir=out_dir,
    slurm_flags="--nodes 1 --gpus-per-node 1",
)


# %%
df_wbm_clean = df_wbm.dropna(subset=MbdKey.init_wyckoff_spglib)

if MbdKey.e_form_dft not in df_wbm_clean:
    raise KeyError(f"{MbdKey.e_form_dft!s} not in {df_wbm_clean.columns=}")
if MbdKey.wyckoff_spglib not in df_wbm_clean:
    raise KeyError(f"{MbdKey.wyckoff_spglib!s} not in {df_wbm_clean.columns=}")


# %%
filters = {
    "created_at": {"$gt": "2022-11-15", "$lt": "2022-11-16"},
    "display_name": {"$regex": "wrenformer-"},
}
runs = wandb.Api().runs(WANDB_PATH, filters=filters)
expected_runs = 10
if len(runs) != expected_runs:
    raise ValueError(
        f"{expected_runs=}, got {len(runs)} filtering {WANDB_PATH=} with {filters=}"
    )

for idx, run in enumerate(runs):
    for key, val in run.config.items():
        if val == runs[0].config[key] or key.startswith(("slurm_", "timestamp")):
            continue
        raise ValueError(
            f"Run configs not identical: runs[{idx}][{key}]={val}, {runs[0][key]=}"
        )

run_params = dict(
    df=dict(shape=str(df_wbm_clean.shape), columns=", ".join(df_wbm_clean)),
    versions={dep: version(dep) for dep in ("aviary", "numpy", "torch")},
    ensemble_size=len(runs),
    task_type=task_type,
    target_col=MbdKey.e_form_dft,
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
data_loader_kwargs = dict(
    input_col=Key.wyckoff,
    target_col=MbdKey.e_form_dft,
    id_col=Key.mat_id,
    embedding_type="wyckoff",
)

data_loader = df_to_in_mem_dataloader(
    df=df_wbm_clean,
    batch_size=1024,
    shuffle=False,  # False is default but best be explicit
    **data_loader_kwargs,
)


# %%
df_pred, ensemble_metrics = predict_from_wandb_checkpoints(
    runs,
    cache_dir=module_dir,
    data_loader=data_loader,
    df=df_wbm_clean,
    model_cls=Wrenformer,
    target_col=MbdKey.e_form_dft,
)
df_pred = df_pred.round(4)

slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")
df_pred.to_csv(f"{out_dir}/{job_name}-preds-{slurm_array_job_id}.csv.gz")


# %%
pred_col = f"{MbdKey.e_form_dft}_pred_ens"
if pred_col not in df_pred:
    raise KeyError(f"{pred_col!s} not in {df_pred.columns=}")
table = wandb.Table(dataframe=df_pred[[MbdKey.e_form_dft, pred_col]].reset_index())


# %%
MAE = ensemble_metrics.MAE.mean()
R2 = ensemble_metrics.R2.mean()

title = f"Wrenformer {task_type} ensemble={len(runs)} {MAE=:.4} {R2=:.4}"

wandb_scatter(table, fields=dict(x=MbdKey.e_form_dft, y=pred_col), title=title)
