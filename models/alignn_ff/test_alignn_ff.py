# %%
from __future__ import annotations

import json
import os
from glob import glob
from importlib.metadata import version

import pandas as pd
import torch
import wandb
from alignn.config import TrainingConfig
from alignn.models.alignn import ALIGNN
from alignn.pretrained import all_models, get_figshare_model
from jarvis.core.graphs import Graph
from pymatgen.core import Structure
from pymatgen.io.jarvis import JarvisAtomsAdaptor
from sklearn.metrics import r2_score
from tqdm import tqdm

from matbench_discovery import Key, Task, today
from matbench_discovery.data import DATA_FILES, df_wbm
from matbench_discovery.plots import wandb_scatter

__author__ = "Philipp Benner, Janosh Riebesell"
__date__ = "2023-07-11"

module_dir = os.path.dirname(__file__)


# %%
n_splits = 100
# model_name = "mp_e_form_alignnn"  # pre-trained by NIST
task_type = Task.IS2RE
device = "cuda" if torch.cuda.is_available() else "cpu"
model_name = f"alignn-ff-wbm-{task_type}"
job_name = f"{model_name}-relaxed-wbm-{task_type}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")
in_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")


if model_name in all_models:  # load pre-trained model
    model = get_figshare_model(model_name)
    pred_col = "e_form_per_atom_alignn_pretrained"
elif os.path.isfile(model_name):
    pred_col = "e_form_per_atom_alignn"
    with open(f"{module_dir}/alignn-config.json") as file:
        config = TrainingConfig(**json.load(file))

    model = ALIGNN(config.model)
    # load trained ALIGNN model
    state_dict = torch.load(model_name, map_location=device)
    model.load_state_dict(state_dict)
    model = model.to(device)
else:
    raise ValueError(
        f"{model_name=} not found, train a model or use pre-trained {list(all_models)}"
    )


# %% Load data
data_path = {
    Task.IS2RE: DATA_FILES.wbm_initial_structures,
    Task.RS2RE: DATA_FILES.wbm_computed_structure_entries,
}[task_type]
# load ALIGNN-FF relaxed structures (TODO fix directory we're loading from)
df_in = pd.concat(map(pd.read_json, glob(f"{module_dir}/data-train-result/*.json.gz")))


# %%
run_params = dict(
    data_path=data_path,
    **{f"{dep}_version": version(dep) for dep in ("alignn", "numpy")},
    model_name=model_name,
    task_type=task_type,
    target_col=Key.e_form,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %% Predict
model.eval()
e_form_preds: dict[str, float] = {}
with torch.no_grad():  # get predictions
    for material_id, structure in tqdm(
        df_in[Key.final_struct].items(),
        total=len(df_in),
        desc=f"Predicting {Key.e_form=} {task_type}",
    ):
        atoms = JarvisAtomsAdaptor.get_atoms(Structure.from_dict(structure))

        atom_graph, line_graph = Graph.atom_dgl_multigraph(atoms)
        e_form = model([atom_graph.to(device), line_graph.to(device)]).item()

        e_form_preds[material_id] = e_form

df_wbm[pred_col] = e_form_preds

df_wbm[pred_col] -= df_wbm.e_correction_per_atom_mp_legacy
df_wbm[pred_col] += df_wbm.e_correction_per_atom_mp2020

if model_name in all_models:
    df_wbm[pred_col].round(4).to_csv(
        f"{module_dir}/{today}-{model_name}-relaxed-wbm-IS2RE.csv.gz"
    )
else:
    df_wbm[pred_col].round(4).to_csv(
        f"{module_dir}/{today}-alignn-relaxed-wbm-IS2RE.csv.gz"
    )


# %%
table = wandb.Table(dataframe=df_wbm[[Key.e_form, pred_col]].reset_index().dropna())

MAE = (df_wbm[Key.e_form] - df_wbm[pred_col]).abs().mean()
R2 = r2_score(df_wbm[Key.e_form], df_wbm[pred_col])
title = f"{model_name} {task_type} {MAE=:.4} {R2=:.4}"
print(title)

wandb_scatter(table, fields=dict(x=Key.e_form, y=pred_col), title=title)
