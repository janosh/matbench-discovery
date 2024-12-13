# %%
import json
import os
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
from pymatviz.enums import Key
from sklearn.metrics import r2_score
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import DataFiles, df_wbm
from matbench_discovery.enums import MbdKey, Task
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell, Philipp Benner"
__date__ = "2023-06-03"

module_dir = os.path.dirname(__file__)


# %%
# model_name = "mp_e_form_alignn"  # pre-trained by NIST (not used for MBD submission)
model_name = DataFiles.alignn_checkpoint.path  # trained by Philipp Benner
task_type = Task.IS2RE
target_col = MbdKey.e_form_dft
input_col = Key.init_struct
device = "cuda" if torch.cuda.is_available() else "cpu"
job_name = f"{model_name}-wbm-{task_type}"
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")


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

slurm_vars = slurm_submit(
    job_name=job_name,
    account="matgen",
    time="11:55:0",
    out_dir=out_dir,
    slurm_flags="--nodes 1 --gpus-per-node 1",
    # pre_cmd is platform specific, remove when running on other systems
    # just left here for reference
    pre_cmd="module load cuda/11.8",
)


# %% Load data
data_path = {
    Task.IS2RE: DataFiles.wbm_initial_structures.path,
    Task.RS2RE: DataFiles.wbm_computed_structure_entries.path,
}[task_type]
input_col = {Task.IS2RE: Key.init_struct, Task.RS2RE: Key.final_struct}[task_type]

df_in = pd.read_json(data_path).set_index(Key.mat_id)

df_in[target_col] = df_wbm[target_col]
if task_type == Task.RS2RE:
    df_in[input_col] = [cse["structure"] for cse in df_in[Key.computed_structure_entry]]
if input_col not in df_in:
    raise KeyError(f"{input_col!s} not in {df_in.columns=}")

df_in[input_col] = [
    JarvisAtomsAdaptor.get_atoms(Structure.from_dict(dct))
    for dct in tqdm(df_in[input_col], leave=False, desc="Converting to JARVIS atoms")
]


# %%
run_params = {
    "data_path": data_path,
    "versions": {dep: version(dep) for dep in ("megnet", "numpy")},
    "model_name": model_name,
    "task_type": task_type,
    "target_col": target_col,
    "df": {"shape": str(df_in.shape), "columns": ", ".join(df_in)},
    "slurm_vars": slurm_vars,
    Key.model_params: sum(  # count trainable params
        p.numel() for p in model.parameters() if p.requires_grad
    ),
}

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %% Predict
model.eval()
e_form_preds: dict[str, float] = {}
with torch.no_grad():  # get predictions
    for material_id, atoms in tqdm(
        df_in[input_col].items(),
        total=len(df_in),
        desc=f"Predicting {target_col=} {task_type}",
    ):
        atom_graph, line_graph = Graph.atom_dgl_multigraph(atoms)
        e_form = model([atom_graph.to(device), line_graph.to(device)]).item()

        e_form_preds[material_id] = e_form

df_wbm[pred_col] = e_form_preds

df_wbm[pred_col] -= df_wbm.e_correction_per_atom_mp_legacy
df_wbm[pred_col] += df_wbm.e_correction_per_atom_mp2020

df_wbm[pred_col].round(4).to_csv(f"{module_dir}/{today}-{model_name}-wbm-IS2RE.csv.gz")


# %%
table = wandb.Table(dataframe=df_wbm[[target_col, pred_col]].reset_index())

MAE = (df_wbm[target_col] - df_wbm[pred_col]).abs().mean()
R2 = r2_score(df_wbm[target_col], df_wbm[pred_col])
title = f"{model_name} {task_type} {MAE=:.4} {R2=:.4}"
print(title)

wandb_scatter(table, fields=dict(x=target_col, y=pred_col), title=title)
