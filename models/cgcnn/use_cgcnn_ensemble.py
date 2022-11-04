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
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
df = pd.read_json(data_path).set_index("material_id", drop=False)
old_len = len(df)
df = df.dropna()  # two missing initial structures
assert len(df) == old_len - 2

df["e_form_per_atom_mp2020_corrected"] = df_wbm.e_form_per_atom_mp2020_corrected

target_col = "e_form_per_atom_mp2020_corrected"
input_col = "initial_structure"
assert target_col in df, f"{target_col=} not in {list(df)}"
assert input_col in df, f"{input_col=} not in {list(df)}"

df[input_col] = [Structure.from_dict(x) for x in tqdm(df[input_col])]

wandb.login()
wandb_api = wandb.Api()
ensemble_id = "cgcnn-e_form-ensemble-1"
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
    df=df,
    target_col=target_col,
    model_class=CrystalGraphConvNet,
    data_loader=data_loader,
)

df.round(6).to_csv(f"{module_dir}/{today}-{ensemble_id}-preds-{target_col}.csv")
