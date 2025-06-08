import json
import os
import pickle as pk

import pandas as pd
import torch
import wandb
from esnet.graphs import (
    PygGraph,
    PygKnowledgeAndStructureDataset,
    prepare_pyg_line_graph_batch,
)
from esnet.models.comformer import iComformer, iComformerConfig
from ignite.handlers import Checkpoint
from pymatgen.core import Structure
from pymatgen.io.jarvis import JarvisAtomsAdaptor
from pymatviz.enums import Key
from sklearn.metrics import r2_score
from torch.utils.data import DataLoader
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import df_wbm
from matbench_discovery.enums import DataFiles, MbdKey, Model, Task
from matbench_discovery.hpc import slurm_submit
from matbench_discovery.plots import wandb_scatter

__author__ = "sl"
__date__ = "2025-06-08"

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


model_name = Model.esnet
task_type = Task.IS2RE
target_col = MbdKey.e_form_dft
job_name = f"{model_name}/{today}-wbm-{task_type}"
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{job_name}")


with open("models/esnet/config.json","rb") as f:
    config = json.load(f)

esnet_config = iComformerConfig(**config["model"])

model = iComformer(esnet_config)
e_form_checkpoint = torch.load(
    "/home/sl/project/ESNet/checkpoint/checkpoint_esnet_259.pt", map_location="cuda:0"
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

atoms = []
for struc in df_in[input_col]:
    lattice = struc.lattice.matrix
    coords = struc.frac_coords
    elements = [site.specie.symbol for site in struc]
    abc = [struc.lattice.a, struc.lattice.b, struc.lattice.c]
    angles = [struc.lattice.alpha, struc.lattice.beta, struc.lattice.gamma]

    atom = {
        "lattice_mat": lattice,
        "coords": coords,
        "elements": elements,
        "abc": abc,
        "angles": angles,
        "cartesian": False,
        "props": [""] * len(elements),
    }

    atoms.append(atom)

df_wbm["atoms"] = atoms

with open("/home/sl/project/ESNet/graphs/RotatE_128_64.pkl", "rb") as f:
    ele2emb = pk.load(f)  # noqa: S301

graphs = []
for struc in df_in[input_col]:
    g = atoms_to_graph(struc)
    graphs.append(g)

graphs = [x.cpu() for x in graphs]


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
data = PygKnowledgeAndStructureDataset(
    df_wbm,
    graphs,
    ele2emb,
    target="e_form_per_atom_mp2020_corrected",
    atom_features="atomic_number",
    line_graph=True,
    id_tag="material_id",
    classification=False,
    neighbor_strategy="k-nearest",
)


wbm_loader = DataLoader(
    data,
    batch_size=64,
    shuffle=False,
    collate_fn=PygKnowledgeAndStructureDataset.collate_line_graph,
    drop_last=False,
    num_workers=1,
)

e_form_values = []

model.eval()
e_form_preds: dict[str, float] = {}
with torch.no_grad():
    for data in wbm_loader:
        x, y = prepare_pyg_line_graph_batch(data, device=device, non_blocking=False)
        output = model(x)
        e_form = output.cpu().numpy().tolist()
        e_form_values.extend(e_form)

pred_col = "e_form_per_atom_esnet"
df_wbm[pred_col] = e_form_values

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
