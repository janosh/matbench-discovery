# %%
import os
from datetime import datetime

import pandas as pd
from aviary.cgcnn.data import CrystalGraphData, collate_batch
from aviary.cgcnn.model import CrystalGraphConvNet
from aviary.core import TaskType
from aviary.train import df_train_test_split, train_model
from pymatgen.core import Structure
from torch.utils.data import DataLoader
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.slurm import slurm_submit

"""
Train a CGCNN ensemble of size n_ens on target_col of data_path.
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-06-13"


# %%
epochs = 300
target_col = "formation_energy_per_atom"
run_name = f"cgcnn-robust-{target_col}"
print(f"{run_name=}")
robust = "robust" in run_name.lower()
n_ens = 10
timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@")[0]
log_dir = f"{os.path.dirname(__file__)}/{today}-{run_name}"

slurm_submit(
    job_name=run_name,
    partition="ampere",
    account="LEE-SL3-GPU",
    time="8:0:0",
    array=f"1-{n_ens}",
    log_dir=log_dir,
    slurm_flags=("--nodes", "1", "--gpus-per-node", "1"),
)


# %%
optimizer = "AdamW"
learning_rate = 3e-4
batch_size = 128
swa_start = None
slurm_array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
task_type: TaskType = "regression"


# %%
data_path = f"{ROOT}/data/mp/2022-08-13-mp-energies.json.gz"
# data_path = f"{ROOT}/data/mp/2022-08-13-mp-energies-1k-samples.json.gz"
print(f"{data_path=}")
df = pd.read_json(data_path).set_index("material_id", drop=False)
df["structure"] = [Structure.from_dict(s) for s in tqdm(df.structure, disable=None)]
assert target_col in df

train_df, test_df = df_train_test_split(df, test_size=0.05)

train_data = CrystalGraphData(train_df, task_dict={target_col: task_type})
train_loader = DataLoader(
    train_data, batch_size=batch_size, shuffle=True, collate_fn=collate_batch
)

test_data = CrystalGraphData(test_df, task_dict={target_col: task_type})
test_loader = DataLoader(
    test_data, batch_size=batch_size, shuffle=False, collate_fn=collate_batch
)

# 1 for regression, n_classes for classification
n_targets = [1 if task_type == "regression" else df[target_col].max() + 1]

model_params = dict(
    n_targets=n_targets,
    elem_emb_len=train_data.elem_emb_len,
    nbr_fea_len=train_data.nbr_fea_dim,
    task_dict={target_col: task_type},  # e.g. {'exfoliation_en': 'regression'}
    robust=robust,
)
model = CrystalGraphConvNet(**model_params)

run_params = dict(
    data_path=data_path,
    batch_size=batch_size,
    train_df=dict(shape=str(train_data.df.shape), columns=", ".join(train_df)),
    test_df=dict(shape=str(test_data.df.shape), columns=", ".join(test_df)),
)


# %%
print(f"\nJob started running {timestamp}")

train_model(
    checkpoint="wandb",  # None | 'local' | 'wandb',
    epochs=epochs,
    learning_rate=learning_rate,
    model_params=model_params,
    model=model,
    optimizer=optimizer,
    run_name=run_name,
    swa_start=swa_start,
    target_col=target_col,
    task_type=task_type,
    test_loader=test_loader,
    timestamp=timestamp,
    train_loader=train_loader,
    wandb_path="janosh/matbench-discovery",
)
