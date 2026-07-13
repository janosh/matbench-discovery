"""Train ALIGNN on Materials Project formation energies."""

# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = [
#   "alignn==2026.5.20",
#   "matbench-discovery",
#   "wandb>=0.27",
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../..", editable = true }
# ///

import json
import os
from importlib.metadata import version

import pandas as pd
import wandb
from alignn.config import TrainingConfig
from alignn.data import get_train_val_loaders
from alignn.train import train_dgl
from pymatgen.core import Structure
from pymatgen.io.jarvis import JarvisAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.enums import DataFiles, Model
from matbench_discovery.hpc import slurm_submit

__author__ = "Philipp Benner, Janosh Riebesell"
__date__ = "2023-06-03"

module_dir = os.path.dirname(__file__)


# %%
model_name = f"{Model.alignn}-mp-e_form"
target_col = Key.formation_energy_per_atom
input_col = "atoms"
job_name = f"{today}-train-{model_name}"


with open(f"{module_dir}/alignn-config.json") as file:
    config = TrainingConfig(**json.load(file))

config.output_dir = out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    account="matgen",
    time="4:0:0",
    out_dir=out_dir,
    slurm_flags="--qos regular --constraint gpu --gpus 1",
)


# %% Load data
df_mp_cse = pd.read_json(DataFiles.mp_computed_structure_entries.path).set_index(
    Key.mat_id
)
df_mp_cse[Key.structure] = [
    Structure.from_dict(cse[Key.structure])
    for cse in tqdm(df_mp_cse.entry, desc="Structures from dict")
]

# load energies
df_in = pd.read_csv(DataFiles.mp_energies.path).set_index(Key.mat_id)
df_in[Key.structure] = df_mp_cse[Key.structure]
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
    alignn_config=config.model_dump(),
)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %%
dataset_array = [
    {
        config.id_tag: material_id,
        input_col: atoms.to_dict(),
        config.target: float(target),
    }
    for material_id, atoms, target in df_in[[input_col, target_col]].itertuples()
]
train_loader, val_loader, test_loader, prepare_batch = get_train_val_loaders(
    dataset=config.dataset,
    dataset_array=dataset_array,
    target=config.target,
    n_train=config.n_train,
    n_val=config.n_val,
    n_test=config.n_test,
    train_ratio=config.train_ratio,
    val_ratio=config.val_ratio,
    test_ratio=config.test_ratio,
    batch_size=config.batch_size,
    atom_features=config.atom_features,
    neighbor_strategy=config.neighbor_strategy,
    standardize=config.atom_features != "cgcnn",
    line_graph=config.compute_line_graph,
    split_seed=config.random_seed,
    id_tag=config.id_tag,
    pin_memory=config.pin_memory,
    workers=config.num_workers,
    save_dataloader=config.save_dataloader,
    use_canonize=config.use_canonize,
    filename=config.filename,
    cutoff=config.cutoff,
    cutoff_extra=config.cutoff_extra,
    max_neighbors=config.max_neighbors,
    three_body_cutoff=config.three_body_cutoff,
    classification_threshold=config.classification_threshold,
    target_multiplication_factor=config.target_multiplication_factor,
    standard_scalar_and_pca=config.standard_scalar_and_pca,
    keep_data_order=config.keep_data_order,
    output_features=config.model.output_features,
    output_dir=config.output_dir,
    use_lmdb=config.use_lmdb,
    read_existing=config.read_existing,
    dtype=config.dtype,
)

# triggers error in alignn/train.py line 1059 in train_dgl()
# f.write("%s, %6f, %6f\n" % (id, target, out_data))
# TypeError: must be real number, not list
config.write_predictions = False

train_hist = train_dgl(
    config,
    train_val_test_loaders=[
        train_loader,
        val_loader,
        test_loader or val_loader,
        prepare_batch,
    ],
)

wandb.log(train_hist)
wandb.save(f"{out_dir}/*")


# %%
df_hist = pd.concat(
    {key: pd.DataFrame(train_hist[key]) for key in ("train", "validation")}
)
table = wandb.Table(dataframe=df_hist)
wandb.log({"table": table})
