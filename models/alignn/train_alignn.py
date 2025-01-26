# %%
import json
import os
from importlib.metadata import version
from typing import Any

import pandas as pd
import torch
import wandb
from alignn.config import TrainingConfig
from alignn.data import StructureDataset, load_graphs
from alignn.train import train_dgl
from pymatgen.core import Structure
from pymatgen.io.jarvis import JarvisAtomsAdaptor
from pymatviz.enums import Key
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import DataFiles
from matbench_discovery.hpc import slurm_submit

__author__ = "Philipp Benner, Janosh Riebesell"
__date__ = "2023-06-03"

module_dir = os.path.dirname(__file__)


# %%
model_name = "alignn-mp-e_form"
target_col = Key.form_energy
input_col = "atoms"
device = "cuda" if torch.cuda.is_available() else "cpu"
job_name = f"train-{model_name}"


pred_col = "e_form_per_atom_alignn"
with open(f"{module_dir}/alignn-config.json") as file:
    config = TrainingConfig(**json.load(file))

config.output_dir = out_dir = os.getenv(
    "SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}"
)

slurm_vars = slurm_submit(
    job_name=job_name,
    account="matgen",
    time="4:0:0",
    out_dir=out_dir,
    slurm_flags="--qos regular --constraint gpu --gpus 1",
)


# %% Load data
df_cse = pd.read_json(DataFiles.mp_computed_structure_entries.path).set_index(
    Key.mat_id
)
df_cse[Key.structure] = [
    Structure.from_dict(cse[Key.structure])
    for cse in tqdm(df_cse.entry, desc="Structures from dict")
]

# load energies
df_in = pd.read_csv(DataFiles.mp_energies.path).set_index(Key.mat_id)
df_in[Key.structure] = df_cse[Key.structure]
if target_col not in df_in:
    raise TypeError(f"{target_col!s} not in {df_in.columns=}")

df_in[input_col] = df_in[Key.structure]
df_in[input_col] = [
    JarvisAtomsAdaptor.get_atoms(struct)
    for struct in tqdm(df_in[Key.structure], desc="Converting to JARVIS atoms")
]


# %%
run_params = dict(
    data_path=DataFiles.mp_energies.path,
    versions={dep: version(dep) for dep in ("alignn", "numpy", "torch", "dgl")},
    model_name=model_name,
    target_col=target_col,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    slurm_vars=slurm_vars,
    alignn_config=config.dict(),
)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %%
df_train, df_val = train_test_split(
    df_in.head(1000).reset_index()[[Key.mat_id, input_col, target_col]],
    test_size=0.05,
    random_state=42,
)


def df_to_loader(
    df: pd.DataFrame,
    *,
    batch_size: int = 128,
    line_graph: bool = True,
    pin_memory: bool = False,
    shuffle: bool = True,
    **kwargs: Any,
) -> DataLoader:
    """Converts a dataframe to a regular PyTorch dataloader for train/val/test.

    Args:
        df (pd.DataFrame): With id, input and target columns
        batch_size (int, optional): Defaults to 128.
        line_graph (bool, optional): Whether to train line (True) or atom (False)  graph
            version of ALIGNN. Defaults to True.
        pin_memory (bool, optional): Whether torch DataLoader should pin memory.
            Defaults to False.
        shuffle (bool, optional): Whether to shuffle the dataset. Defaults to True.
        **kwargs: Additional arguments to pass to the StructureDataset

    Returns:
        DataLoader: PyTorch data loader
    """
    graphs = load_graphs(
        df, neighbor_strategy=config.neighbor_strategy, use_canonize=config.use_canonize
    )
    dataset = StructureDataset(
        df.reset_index(drop=True),
        graphs,
        target=target_col,
        line_graph=line_graph,
        atom_features=config.atom_features,
        id_tag=Key.mat_id,
        **kwargs,
    )
    collate_fn = getattr(dataset, f"collate{'_line' if line_graph else ''}_graph")

    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        collate_fn=collate_fn,
        pin_memory=pin_memory,
    )


train_loader, val_loader = df_to_loader(df_train), df_to_loader(df_val, shuffle=False)


# %%
prepare_batch = train_loader.dataset.prepare_batch
# triggers error in alignn/train.py line 1059 in train_dgl()
# f.write("%s, %6f, %6f\n" % (id, target, out_data))
# TypeError: must be real number, not list
config.write_predictions = False

train_hist = train_dgl(
    config,
    train_val_test_loaders=[train_loader, val_loader, val_loader, prepare_batch],
)

wandb.log(train_hist)
wandb.save(f"{out_dir}/*")


# %%
df_hist = pd.concat(
    {key: pd.DataFrame(train_hist[key]) for key in ("train", "validation")}
)
table = wandb.Table(dataframe=df_hist)
wandb.log({"table": table})
