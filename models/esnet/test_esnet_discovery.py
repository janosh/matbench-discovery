import os
import pickle as pk
from collections import Counter
from typing import Any

import numpy as np
import pandas as pd
import torch
import wandb
from esnet.graphs import PygGraph, PygStructureDataset
from esnet.load_triples import Triples
from esnet.models.comformer import iComformer
from ignite.handlers import Checkpoint
from pymatgen.core import Structure
from pymatgen.io.jarvis import JarvisAtomsAdaptor
from pymatviz.enums import Key
from sklearn.metrics import r2_score
from torch_geometric.data import Data
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import df_wbm
from matbench_discovery.enums import DataFiles, MbdKey, Model, Task
from matbench_discovery.hpc import slurm_submit
from matbench_discovery.plots import wandb_scatter

__author__ = "sl"
__date__ = "2025-05-13"

device = "cuda" if torch.cuda.is_available() else "cpu"
adaptor = JarvisAtomsAdaptor()


def atoms_to_graph(atoms: Structure) -> PygGraph:
    """Convert structure to Atom."""
    structure = adaptor.get_atoms(atoms)
    return PygGraph.atom_dgl_multigraph(
        structure,
        neighbor_strategy="k-nearest",
        cutoff=8.0,
        atom_features="atomic_number",
        max_neighbors=25,
        compute_line_graph=False,
        use_canonize=True,
        use_lattice=True,
        use_angle=False,
    )


def graph_build(atoms: Structure) -> Data:
    g = atoms_to_graph(atoms)
    features = PygStructureDataset.get_attribute_lookup("cgcnn")
    z = g.x
    g.atomic_number = z
    z = z.type(torch.IntTensor).squeeze()
    f = torch.tensor(features[z]).type(torch.FloatTensor)
    if g.x.size(0) == 1:
        f = f.unsqueeze(0)
    g.x = f
    g.batch = torch.zeros(1, dtype=torch.int64)
    return g


with open("/home/sl/project/ESNet/graphs/RotatE_128_64.pkl", "rb") as f:
    rotate_emb = pk.load(f)  # noqa: S301


def load_elements_knowledge(elements: list[str], rotate_emb: dict[str, Any]) -> Data:
    data = Triples()

    entity = rotate_emb["entity"]
    element_counts = Counter(elements)
    total_elements = len(elements)
    element_ratios = {
        element: count / total_elements for element, count in element_counts.items()
    }
    template = entity[0]
    features = np.zeros((len(element_counts), len(np.array(template))))

    attributes = []

    for z, (element, ratios) in enumerate(element_ratios.items()):
        if element in data.entities:
            attribute_id = [h for (r, h) in data.t2rh[data.entity2id[element]]]
            attributes.extend(attribute_id)

        key = data.entity2id[element]
        x = entity[int(key)]
        if x is not None:
            features[z, :] = x * ratios

    attributes = sorted(set(attributes))
    attr_features = np.zeros((len(attributes), len(np.array(template))))

    for t, key in enumerate(attributes):
        attr_features[t, :] = entity[int(key)]

    f = torch.tensor(features).type(torch.get_default_dtype())
    attr_f = torch.tensor(attr_features).type(torch.get_default_dtype())
    node_emb = torch.cat((f, attr_f), 0)
    knowledge_features = Data(x=node_emb)
    knowledge_features.batch = torch.zeros(1, dtype=torch.int64)

    return knowledge_features


model_name = Model.esnet
task_type = Task.IS2RE
target_col = MbdKey.e_form_dft
job_name = f"{model_name}/{today}-wbm-{task_type}"
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{job_name}")

model = iComformer()
e_form_checkpoint = torch.load(
    "/home/sl/project/ESNet/checkpoint/checkpoint_eform_500.pt", map_location="cuda:0"
)
to_load = {
    "model": model,
}
Checkpoint.load_objects(to_load=to_load, checkpoint=e_form_checkpoint)
model = model.to(device)


slurm_vars = slurm_submit(
    job_name=job_name,
    account="matgen",
    time="12:0:0",
    out_dir=out_dir,
    slurm_flags="--nodes 1 --gpus-per-node 1",
    pre_cmd="module load cuda/12.1",
)


# %% Load data
data_path = {
    Task.IS2RE: DataFiles.wbm_initial_structures.path,
    Task.RS2RE: DataFiles.wbm_computed_structure_entries.path,
}[task_type]
input_col = {Task.IS2RE: "initial_structure", Task.RS2RE: Key.final_struct}[task_type]

df_in = pd.read_json(data_path, lines=True).set_index(Key.mat_id)

if input_col not in df_in:
    raise TypeError(f"{input_col!s} not in {df_in.columns=}")

df_in[MbdKey.e_form_dft] = df_wbm[MbdKey.e_form_dft]
if task_type == Task.RS2RE:
    df_in[input_col] = [cse["structure"] for cse in df_in[Key.computed_structure_entry]]

df_in[input_col] = [
    Structure.from_dict(dct) for dct in tqdm(df_in[input_col], disable=None)
]


# %%
run_params = {
    "data_path": data_path,
    # "versions": {dep: version(dep) for dep in ("megnet", "numpy")},
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
# eform 均值与方差
std_eform = 1.073214
mean_eform = -1.652392

model.eval()
e_form_preds: dict[str, float] = {}
with torch.no_grad():  # get predictions
    for material_id, atoms in tqdm(
        df_in[input_col].items(),
        total=len(df_in),
        desc=f"Predicting {target_col=} {task_type}",
    ):
        elements = [site.specie.symbol for site in atoms.sites]
        g = graph_build(atoms)
        kge = load_elements_knowledge(elements, rotate_emb)
        out_data = model([g.to(device), g.to(device), kge.to(device), g.to(device)])
        out_data = out_data.cpu().numpy()
        e_form = out_data * std_eform + mean_eform
        e_form_preds[material_id] = e_form

pred_col = "e_form_per_atom_esnet"
df_wbm[pred_col] = e_form_preds

df_wbm[pred_col] -= df_wbm.e_correction_per_atom_mp_legacy
df_wbm[pred_col] += df_wbm.e_correction_per_atom_mp2020
df_wbm[pred_col].round(4).to_csv(f"{module_dir}/{today}-{model_name}-wbm-IS2RE.csv.gz")


# %%
table = wandb.Table(dataframe=df_wbm[[target_col, pred_col]].reset_index())

MAE = (df_wbm[target_col] - df_wbm[pred_col]).abs().mean()
R2 = r2_score(df_wbm[target_col], df_wbm[pred_col])

print("MAE:", MAE)
print("R2:", R2)

title = f"{model_name} {task_type} {MAE=:.4} {R2=:.4}"
print("title:", title)


wandb_scatter(table, fields=dict(x=target_col, y=pred_col), title=title)
