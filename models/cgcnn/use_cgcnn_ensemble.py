# %%
from __future__ import annotations

import os
from datetime import datetime

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
ensemble_id = "cgcnn-e_form-ensemble-1"
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
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
df = pd.read_json(data_path).set_index("material_id", drop=False)
old_len = len(df)
no_init_structs = df.query("initial_structure.isnull()").index
df = df.dropna()  # two missing initial structures
assert len(df) == old_len - 2

assert all(df.index == df_wbm.drop(index=no_init_structs).index)

target_col = "e_form_per_atom_mp2020_corrected"
df[target_col] = df_wbm[target_col]
input_col = "initial_structure"
assert target_col in df, f"{target_col=} not in {list(df)}"
assert input_col in df, f"{input_col=} not in {list(df)}"

df[input_col] = [Structure.from_dict(x) for x in tqdm(df[input_col], disable=None)]

wandb.login()
wandb_api = wandb.Api()
runs = wandb_api.runs(
    "janosh/matbench-discovery", filters={"tags": {"$in": [ensemble_id]}}
)

assert len(runs) == 10, f"Expected 10 runs, got {len(runs)} for {ensemble_id=}"

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

df.round(6).to_csv(f"{module_dir}/{today}-{run_name}-preds.csv", index=False)
