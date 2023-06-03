# %%
from __future__ import annotations

import os

import pandas as pd
import torch
from alignn.pretrained import get_figshare_model
from jarvis.core.graphs import Graph
from pymatgen.core import Structure
from pymatgen.io.jarvis import JarvisAtomsAdaptor
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import DATA_FILES, df_wbm

__author__ = "Janosh Riebesell"
__date__ = "2023-06-03"

module_dir = os.path.dirname(__file__)


# %%
model = get_figshare_model("mp_e_form_alignnn")
task_type = "IS2RE"
target_col = "e_form_per_atom_mp2020_corrected"
input_col = "initial_structure"
id_col = "material_id"
pred_col = "e_form_per_atom_alignn"
device = "cuda" if torch.cuda.is_available() else "cpu"


# %% Load data
data_path = {
    "IS2RE": DATA_FILES.wbm_initial_structures,
    "RS2RE": DATA_FILES.wbm_computed_structure_entries,
}[task_type]
input_col = {"IS2RE": "initial_structure", "RS2RE": "relaxed_structure"}[task_type]

df = pd.read_json(data_path).set_index(id_col)

df[target_col] = df_wbm[target_col]
if task_type == "RS2RE":
    df[input_col] = [x["structure"] for x in df.computed_structure_entry]
assert input_col in df, f"{input_col=} not in {list(df)}"

df[input_col] = [
    JarvisAtomsAdaptor.get_atoms(Structure.from_dict(x))
    for x in tqdm(df[input_col], disable=None)
]


# %% Compute test result
model.eval()
e_form_preds: dict[str, float] = {}
with torch.no_grad():  # get predictions
    for material_id, atoms in tqdm(df[input_col].items()):
        atom_graph, line_graph = Graph.atom_dgl_multigraph(atoms)
        e_form = model([atom_graph.to(device), line_graph.to(device)]).item()

        e_form_preds[material_id] = e_form

df_wbm[pred_col] = e_form_preds

df_wbm[pred_col].round(4).to_csv(f"{module_dir}/{today}-alignn-wbm-IS2RE.csv")
