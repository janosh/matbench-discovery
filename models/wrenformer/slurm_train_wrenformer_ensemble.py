# %%
import os
from datetime import datetime

import pandas as pd
from aviary.train import df_train_test_split, train_wrenformer

from matbench_discovery import ROOT
from matbench_discovery.slurm import slurm_submit_python

"""
Train a Wrenformer ensemble of size n_folds on target_col of data_path.
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-08-13"


# %%
epochs = 300
data_path = f"{ROOT}/data/mp/2022-08-13-mp-energies.json.gz"
target_col = "formation_energy_per_atom"
# data_path = f"{ROOT}/data/2022-08-25-m3gnet-trainset-mp-2021-struct-energy.json.gz"
# target_col = "mp_energy_per_atom"
data_name = "m3gnet-trainset" if "m3gnet" in data_path else "mp"
run_name = f"wrenformer-robust-{data_name}-{target_col}-{epochs=}"
n_folds = 10
timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@")[0]
dataset = "mp"
log_dir = f"{os.path.dirname(__file__)}/{dataset}/{today}-{run_name}"

slurm_submit_python(
    job_name=run_name,
    partition="ampere",
    account="LEE-SL3-GPU",
    time="8:0:0",
    array=f"1-{n_folds}",
    log_dir=log_dir,
    slurm_flags=("--nodes", "1", "--gpus-per-node", "1"),
)


# %%
learning_rate = 3e-4
batch_size = 128
slurm_array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
input_col = "wyckoff_spglib"

print(f"Job started running {timestamp}")
print(f"{run_name=}")
print(f"{data_path=}")

df = pd.read_json(data_path).set_index("material_id", drop=False)
assert target_col in df, f"{target_col=} not in {list(df)}"
assert input_col in df, f"{input_col=} not in {list(df)}"
train_df, test_df = df_train_test_split(df, test_size=0.05)

run_params = dict(
    data_path=data_path,
    batch_size=batch_size,
    train_df=dict(shape=train_df.shape, columns=", ".join(train_df)),
    test_df=dict(shape=test_df.shape, columns=", ".join(test_df)),
)

train_wrenformer(
    run_name=run_name,
    train_df=train_df,
    test_df=test_df,
    target_col=target_col,
    task_type="regression",
    timestamp=timestamp,
    # folds=(n_folds, slurm_array_task_id),
    epochs=epochs,
    checkpoint="wandb",  # None | 'local' | 'wandb',
    input_col=input_col,
    learning_rate=learning_rate,
    batch_size=batch_size,
    wandb_path="janosh/matbench-discovery",
    run_params=run_params,
)
