# %%
from __future__ import annotations

import os
from datetime import datetime
from importlib.metadata import version

import pandas as pd
import wandb
from aviary.cgcnn.data import CrystalGraphData, collate_batch
from aviary.cgcnn.model import CrystalGraphConvNet
from aviary.deploy import predict_from_wandb_checkpoints
from pymatgen.core import Structure
from torch.utils.data import DataLoader
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.plot_scripts import df_wbm
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

"""
Script that downloads checkpoints for an ensemble of Wrenformer models trained on the MP
formation energies, then makes predictions on some dataset, prints ensemble metrics and
stores predictions to CSV.
"""

today = f"{datetime.now():%Y-%m-%d}"
log_dir = f"{os.path.dirname(__file__)}/{today}-test"
job_name = "test-cgcnn-ensemble"

slurm_vars = slurm_submit(
    job_name=job_name,
    partition="ampere",
    account="LEE-SL3-GPU",
    time=(slurm_max_job_time := "2:0:0"),
    log_dir=log_dir,
    slurm_flags=("--nodes", "1", "--gpus-per-node", "1"),
)


# %%
task_type = "IS2RE"
if task_type == "IS2RE":
    data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
elif task_type == "RS2RE":
    data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-cses.json.bz2"
df = pd.read_json(data_path).set_index("material_id", drop=False)

target_col = "e_form_per_atom_mp2020_corrected"
df[target_col] = df_wbm[target_col]
input_col = "initial_structure"
assert target_col in df, f"{target_col=} not in {list(df)}"
assert input_col in df, f"{input_col=} not in {list(df)}"

df[input_col] = [Structure.from_dict(x) for x in tqdm(df[input_col], disable=None)]

filters = {
    "$and": [{"created_at": {"$gt": "2022-11-22", "$lt": "2022-11-23"}}],
    "display_name": {"$regex": "^cgcnn-robust"},
}
wandb.login()
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
    slurm_vars=slurm_vars | dict(slurm_max_job_time=slurm_max_job_time),
)

slurm_job_id = os.environ.get("SLURM_JOB_ID", "debug")
wandb.init(
    project="matbench-discovery", name=f"{job_name}-{slurm_job_id}", config=run_params
)

cg_data = CrystalGraphData(
    df, task_dict={target_col: "regression"}, structure_col=input_col
)
data_loader = DataLoader(
    cg_data, batch_size=1024, shuffle=False, collate_fn=collate_batch
)
df, ensemble_metrics = predict_from_wandb_checkpoints(
    runs,
    # dropping isolated-atom structs means len(cg_data.df) < len(df)
    df=cg_data.df.reset_index(drop=True).drop(columns=input_col),
    target_col=target_col,
    model_cls=CrystalGraphConvNet,
    data_loader=data_loader,
)

df.to_csv(f"{log_dir}/{today}-{job_name}-preds.csv", index=False)
table = wandb.Table(dataframe=df)


# %%
pred_col = f"{target_col}_pred_ens"
MAE = ensemble_metrics["MAE"]
R2 = ensemble_metrics["R2"]

title = rf"CGCNN {task_type} ensemble={len(runs)} {MAE=:.4} {R2=:.4}"
print(title)

scatter_plot = wandb.plot_table(
    vega_spec_name="janosh/scatter-parity",
    data_table=table,
    fields=dict(x=target_col, y=pred_col, title=title),
)

wandb.log({"true_pred_scatter": scatter_plot})
