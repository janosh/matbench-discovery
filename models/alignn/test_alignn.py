# %%
from __future__ import annotations

import csv
import os
from typing import Any

import torch
from alignn.config import TrainingConfig
from alignn.models.alignn import ALIGNN
from jarvis.core.atoms import Atoms
from jarvis.core.graphs import Graph
from jarvis.db.jsonutils import loadjson
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import df_wbm

__author__ = "Philipp Benner"
__date__ = "2023-06-02"

module_dir = os.path.dirname(__file__)
pred_col = "e_form_per_atom_alignn"


# %% ALIGNN config
device = "cuda" if torch.cuda.is_available() else "cpu"
config = loadjson(f"{module_dir}/alignn-config.json")
config = TrainingConfig(**config)


# %% Load test data
def load_data_directory(basename: str) -> list[dict[str, Any]]:
    """Load ASE atoms, material IDs and target values from a directory."""
    id_prop_dat = os.path.join(basename, "id_prop.csv")

    with open(id_prop_dat) as file:
        reader = csv.reader(file)
        data = list(reader)

    dataset = []
    for sample in data:
        filename = os.path.join(basename, sample[0])

        info = {}
        info["atoms"] = Atoms.from_poscar(filename)
        info["jid"] = sample[0]
        info["target"] = float(sample[1])

        dataset.append(info)

    return dataset


# %% Compute test result
model = ALIGNN(config.model)
# load trained ALIGNN model
state_dict = torch.load(
    f"{module_dir}/data-train-result/best-model.pth", map_location=device
)
model.load_state_dict(state_dict)
model = model.to(device)

dataset = load_data_directory("data-test-wbm")

model.eval()
predictions = []
with torch.no_grad():  # get predictions
    for datum in tqdm(dataset):
        atom_graph, line_graph = Graph.atom_dgl_multigraph(
            datum["atoms"],
            cutoff=config.cutoff,
            atom_features=config.atom_features,
            max_neighbors=config.max_neighbors,
        )
        y_hat = model([atom_graph.to(device), line_graph.to(device)]).item()

        predictions.append(y_hat)

df_wbm[pred_col] = predictions

df_wbm[pred_col].round(4).to_csv(f"{module_dir}/{today}-alignn-wbm-IS2RE.csv")
