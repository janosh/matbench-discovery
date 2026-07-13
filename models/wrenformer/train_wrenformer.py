"""Train a Wrenformer ensemble on Materials Project formation energies."""

# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = [
#   "aviary-models==1.2.1",
#   "matbench-discovery",
#   "torch",
#   "wandb>=0.27",
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../..", editable = true }
# ///

import os
import shutil
from importlib.metadata import version

import pandas as pd
import torch
import wandb
from aviary import ROOT as AVIARY_ROOT
from aviary.utils import train_ensemble
from aviary.wrenformer.data import df_to_in_mem_dataloader
from aviary.wrenformer.model import Wrenformer
from pymatviz.enums import Key

from matbench_discovery import timestamp, today
from matbench_discovery.enums import DataFiles, MbdKey, Model
from matbench_discovery.hpc import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-13"


# %%
epochs = 300
data_path = DataFiles.mp_energies.path
target_col = Key.formation_energy_per_atom
# data_path = f"{ROOT}/data/2022-08-25-m3gnet-trainset-mp-2021-struct-energy.json.gz"
# target_col = "mp_energy_per_atom"
data_name = "m3gnet-trainset" if "m3gnet" in data_path else "mp"
ensemble_size = 10
job_name = f"{today}-train-{Model.wrenformer}-ens={ensemble_size}-robust-{data_name}"
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{job_name}")


slurm_vars = slurm_submit(
    job_name=job_name,
    account="matgen",
    time="8:0:0",
    array=f"1-{ensemble_size}",
    out_dir=out_dir,
    slurm_flags="--nodes 1 --gpus-per-node 1",
)


# %%
learning_rate = 3e-4
batch_size = 128
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "1"))
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print(f"\nJob {job_name} started {timestamp}")

df_data = pd.read_csv(data_path, keep_default_na=False, na_values=[]).set_index(
    Key.mat_id, drop=False
)

if target_col not in df_data:
    raise TypeError(f"{target_col!s} not in {df_data.columns=}")
if MbdKey.protostructure_spglib not in df_data:
    raise TypeError(f"{MbdKey.protostructure_spglib!s} not in {df_data.columns=}")

# Match the deterministic split used by aviary v0.1.0's df_train_test_split().
df_shuffled = df_data.sample(frac=1, random_state=0)
train_df = df_shuffled.sample(frac=0.95, random_state=0)
test_df = df_shuffled.drop(train_df.index)

loss_dict = {target_col: "L1"}
data_loader_kwargs = dict(
    input_col=MbdKey.protostructure_spglib,
    target_col=target_col,
    id_col=Key.mat_id,
    embedding_type="protostructure",
    device=device,
)
train_loader = df_to_in_mem_dataloader(
    train_df,
    batch_size=batch_size,
    shuffle=True,
    **data_loader_kwargs,
)
test_loader = df_to_in_mem_dataloader(
    test_df,
    batch_size=batch_size * 16,
    shuffle=False,
    **data_loader_kwargs,
)

model_name = f"{Model.wrenformer}-robust-{data_name}"
model_params = dict(
    task_dict={target_col: "regression"},
    robust=True,
    n_targets=[1],
    n_features=train_loader.tensors[0][0].shape[-1],
)
setup_params = dict(
    optim="AdamW",
    learning_rate=learning_rate,
    weight_decay=0.01,
    momentum=0.9,
    device=device,
)
restart_params = dict(resume=None, fine_tune=None, transfer=None)
run_params = dict(
    data_path=data_path,
    epochs=epochs,
    versions={
        dep: version(dep) for dep in ("aviary-models", "numpy", "torch", "wandb")
    },
    batch_size=batch_size,
    train_df=dict(shape=train_df.shape, columns=", ".join(train_df)),
    test_df=dict(shape=test_df.shape, columns=", ".join(test_df)),
    slurm_vars=slurm_vars,
    task_type="regression",
    loss_dict=loss_dict,
    model_params=model_params,
    setup_params={**setup_params, "device": str(device)},
    restart_params=restart_params,
)


# %%
# Pre-initializing keeps Aviary's W&B logging in the public benchmark project.
wandb.init(
    entity="janosh",
    project="matbench-discovery",
    name=f"{job_name}-{slurm_array_task_id}",
    config=run_params,
)
train_ensemble(
    model_class=Wrenformer,
    model_name=model_name,
    # folds=(ensemble_size, slurm_array_task_id),
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

# The legacy script used checkpoint="wandb" (alternatives: None or "local").
aviary_checkpoint = (
    f"{AVIARY_ROOT}/models/{model_name}/checkpoint-r{slurm_array_task_id}.pth.tar"
)
os.makedirs(out_dir, exist_ok=True)
checkpoint_path = f"{out_dir}/checkpoint-r{slurm_array_task_id}.pth.tar"
shutil.copy2(aviary_checkpoint, checkpoint_path)
wandb.save(checkpoint_path, base_path=out_dir)
wandb.finish()
