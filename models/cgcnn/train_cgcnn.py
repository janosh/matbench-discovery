"""Train a CGCNN ensemble on Materials Project formation energies."""

# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = [
#   "aviary-models==1.2.1",
#   "matbench-discovery",
#   "numpy>=2,<3",
#   "torch",
#   "tqdm>=4.67.3",
#   "wandb>=0.27",
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../..", editable = true }
# ///

import os
import shutil
from importlib.metadata import version

import numpy as np
import pandas as pd
import torch
import wandb
from aviary import ROOT as AVIARY_ROOT
from aviary.cgcnn.data import CrystalGraphData, collate_batch
from aviary.cgcnn.model import CrystalGraphConvNet
from aviary.utils import train_ensemble
from pymatgen.core import Structure
from pymatviz.enums import Key
from torch.utils.data import DataLoader
from tqdm import tqdm, trange

from matbench_discovery import timestamp, today
from matbench_discovery.enums import DataFiles
from matbench_discovery.hpc import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-06-13"

np_rng = np.random.default_rng(seed=0)  # ensure reproducible structure perturbations


def perturb_structure(struct: Structure, gamma: float = 1.5) -> Structure:
    """Perturb atomic coordinates for CGCNN+P training set augmentation."""
    perturbed = struct.copy()
    for site in perturbed:
        magnitude = np_rng.weibull(gamma)
        vector = np_rng.normal(size=3)
        norm = np.linalg.norm(vector)
        vector = vector / norm if norm > np.finfo(float).eps else np.array([1.0, 0, 0])
        site.coords += vector * magnitude
        site.to_unit_cell(in_place=True)

    return perturbed


# %%
epochs = 300
target_col = Key.formation_energy_per_atom
input_col = Key.structure
# 0 for no perturbation, n>1 means train on n perturbations of each crystal
# in the training set all assigned the same original target energy
n_perturb = 0
model_name = f"cgcnn-robust-{n_perturb}"
job_name = f"{today}-train-{model_name}"
print(f"{job_name=}")
ensemble_size = 10
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    account="matgen",
    time="11:55:0",
    array=f"1-{ensemble_size}",
    out_dir=out_dir,
    slurm_flags="--nodes 1 --gpus-per-node 1",
)


# %%
learning_rate = 3e-4
batch_size = 128
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "1"))
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# %%
data_path = DataFiles.mp_energies.path
df_in = pd.read_csv(data_path).set_index(Key.mat_id)

df_mp_cse = pd.read_json(DataFiles.mp_computed_structure_entries.path).set_index(
    Key.mat_id
)
df_in[input_col] = [
    Structure.from_dict(computed_entry[input_col])
    for computed_entry in tqdm(df_mp_cse.entry, desc="Structures from dict")
]

if target_col not in df_in:
    raise TypeError(f"{target_col!s} not in {df_in.columns=}")

df_aug = df_in.copy()
structures = df_aug.pop(input_col)
material_ids = df_aug.index
for perturb_idx in trange(n_perturb, desc="Generating perturbed structures"):
    df_aug[input_col] = [perturb_structure(struct) for struct in structures]
    df_aug.index = [
        f"{material_id}-aug={perturb_idx + 1}" for material_id in material_ids
    ]
    df_in = pd.concat([df_in, df_aug])

del df_aug

# Match the deterministic split used by aviary v0.1.0's df_train_test_split().
df_shuffled = df_in.sample(frac=1, random_state=0)
train_df = df_shuffled.sample(frac=0.95, random_state=0)
test_df = df_shuffled.drop(train_df.index)
task_dict = {target_col: "regression"}

print(f"{train_df.shape=}")
train_data = CrystalGraphData(train_df, task_dict=task_dict)
train_loader = DataLoader(
    train_data, batch_size=batch_size, shuffle=True, collate_fn=collate_batch
)

print(f"{test_df.shape=}")
test_data = CrystalGraphData(test_df, task_dict=task_dict)
test_loader = DataLoader(
    test_data, batch_size=batch_size, shuffle=False, collate_fn=collate_batch
)

model_params = dict(
    n_targets=train_data.n_targets,
    task_dict=task_dict,
    robust=True,
)
setup_params = dict(
    optim="AdamW",
    learning_rate=learning_rate,
    weight_decay=0.01,
    momentum=0.9,
    device=device,
)
restart_params = dict(resume=None, fine_tune=None, transfer=None)
loss_dict = {target_col: "L1"}
run_name = f"{job_name}-{slurm_array_task_id}"
run_params = dict(
    data_path=data_path,
    epochs=epochs,
    batch_size=batch_size,
    versions={
        dep: version(dep) for dep in ("aviary-models", "numpy", "torch", "wandb")
    },
    train_df=dict(shape=str(train_data.df.shape), columns=", ".join(train_df)),
    test_df=dict(shape=str(test_data.df.shape), columns=", ".join(test_df)),
    slurm_vars=slurm_vars,
    n_perturb=n_perturb,
    input_col=input_col,
    task_type="regression",
    loss_dict=loss_dict,
    model_params=model_params,
    setup_params={**setup_params, "device": str(device)},
    restart_params=restart_params,
)


# %%
print(f"\nJob {job_name} started {timestamp}")

# Pre-initializing keeps Aviary's W&B logging in the public benchmark project.
wandb.init(
    entity="janosh", project="matbench-discovery", name=run_name, config=run_params
)
train_ensemble(
    model_class=CrystalGraphConvNet,
    model_name=model_name,
    run_id=slurm_array_task_id,
    ensemble_folds=1,
    epochs=epochs,
    patience=None,
    train_loader=train_loader,
    val_loader=test_loader,
    log="wandb",
    setup_params=setup_params,
    restart_params=restart_params,
    model_params=model_params,
    loss_dict=loss_dict,
)

aviary_checkpoint = (
    f"{AVIARY_ROOT}/models/{model_name}/checkpoint-r{slurm_array_task_id}.pth.tar"
)
os.makedirs(out_dir, exist_ok=True)
checkpoint_path = f"{out_dir}/checkpoint-r{slurm_array_task_id}.pth.tar"
shutil.copy2(aviary_checkpoint, checkpoint_path)
wandb.save(checkpoint_path, base_path=out_dir)
wandb.finish()
