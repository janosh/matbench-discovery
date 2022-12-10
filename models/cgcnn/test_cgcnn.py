# %%
from __future__ import annotations

import os
import sys
from importlib.metadata import version

import pandas as pd
import wandb
from aviary.cgcnn.data import CrystalGraphData, collate_batch
from aviary.cgcnn.model import CrystalGraphConvNet
from aviary.deploy import predict_from_wandb_checkpoints
from pymatgen.core import Structure
from torch.utils.data import DataLoader
from tqdm import tqdm

from matbench_discovery import CHECKPOINT_DIR, DEBUG, ROOT, WANDB_PATH, today
from matbench_discovery.data import df_wbm
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"

"""
Download WandB checkpoints for an ensemble of CGCNN models trained on all MP
formation energies, then makes predictions on some dataset, prints ensemble metrics and
saves predictions to CSV.
"""

task_type = "IS2RE"
debug = "slurm-submit" in sys.argv
job_name = f"test-cgcnn-wbm-{task_type}{'-debug' if DEBUG else ''}"
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
if task_type == "IS2RE":
    data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
    input_col = "initial_structure"
elif task_type == "RS2RE":
    data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-computed-structure-entries.json.bz2"
    input_col = "relaxed_structure"
else:
    raise ValueError(f"Unexpected {task_type=}")

df = pd.read_json(data_path).set_index("material_id")

target_col = "e_form_per_atom_mp2020_corrected"
df[target_col] = df_wbm[target_col]
if task_type == "RS2RE":
    df[input_col] = [x["structure"] for x in df.computed_structure_entry]
assert input_col in df, f"{input_col=} not in {list(df)}"

df[input_col] = [Structure.from_dict(x) for x in tqdm(df[input_col], disable=None)]

filters = {
    "created_at": {"$gt": "2022-12-03", "$lt": "2022-12-04"},
    "display_name": {"$regex": "^train-cgcnn-robust-augment=3-"},
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
    ensemble_size=len(runs),
    task_type=task_type,
    target_col=target_col,
    input_col=input_col,
    wandb_run_filters=filters,
    slurm_vars=slurm_vars,
)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)

cg_data = CrystalGraphData(
    df, task_dict={target_col: "regression"}, structure_col=input_col
)
data_loader = DataLoader(
    cg_data, batch_size=1024, shuffle=False, collate_fn=collate_batch
)
df_preds, ensemble_metrics = predict_from_wandb_checkpoints(
    runs,
    # dropping isolated-atom structs means len(cg_data.df) < len(df)
    cache_dir=CHECKPOINT_DIR,
    df=cg_data.df.drop(columns=input_col),
    target_col=target_col,
    model_cls=CrystalGraphConvNet,
    data_loader=data_loader,
)

slurm_job_id = os.environ.get("SLURM_JOB_ID", "debug")
df_preds.round(4).to_csv(f"{out_dir}/{job_name}-preds-{slurm_job_id}.csv", index=False)
pred_col = f"{target_col}_pred_ens"
assert pred_col in df, f"{pred_col=} not in {list(df)}"
table = wandb.Table(dataframe=df_preds[[target_col, pred_col]].reset_index())


# %%
MAE = ensemble_metrics.MAE.mean()
R2 = ensemble_metrics.R2.mean()

title = f"CGCNN {task_type} ensemble={len(runs)} {MAE=:.4} {R2=:.4}"

wandb_scatter(table, fields=dict(x=target_col, y=pred_col), title=title)
