"""Train a Wrenformer ensemble on target_col of data_path."""


# %%
import os
from importlib.metadata import version

import pandas as pd
from aviary.train import df_train_test_split, train_wrenformer

from matbench_discovery import DEBUG, WANDB_PATH, timestamp, today
from matbench_discovery.data import DATA_FILES
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-13"


# %%
epochs = 300
data_path = DATA_FILES.mp_energies
target_col = "formation_energy_per_atom"
# data_path = f"{ROOT}/data/2022-08-25-m3gnet-trainset-mp-2021-struct-energy.json.gz"
# target_col = "mp_energy_per_atom"
data_name = "m3gnet-trainset" if "m3gnet" in data_path else "mp"
job_name = f"train-wrenformer-robust-{data_name}{'-debug' if DEBUG else ''}"
ensemble_size = 10
dataset = "mp"
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")


slurm_vars = slurm_submit(
    job_name=job_name,
    partition="ampere",
    account="LEE-SL3-GPU",
    time="8:0:0",
    array=f"1-{ensemble_size}",
    out_dir=out_dir,
    slurm_flags="--nodes 1 --gpus-per-node 1",
)


# %%
learning_rate = 3e-4
batch_size = 128
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
input_col = "wyckoff_spglib"

print(f"\nJob started running {timestamp}")
print(f"{job_name=}")
print(f"{data_path=}")

df = pd.read_csv(data_path).set_index("material_id", drop=False)

assert target_col in df, f"{target_col=} not in {list(df)}"
assert input_col in df, f"{input_col=} not in {list(df)}"
train_df, test_df = df_train_test_split(df, test_size=0.05)

run_params = dict(
    data_path=data_path,
    **{f"{dep}_version": version(dep) for dep in ("aviary", "numpy", "torch")},
    batch_size=batch_size,
    train_df=dict(shape=train_df.shape, columns=", ".join(train_df)),
    test_df=dict(shape=test_df.shape, columns=", ".join(test_df)),
    slurm_vars=slurm_vars,
)

train_wrenformer(
    run_name=f"{job_name}-{slurm_array_task_id}",
    train_df=train_df,
    test_df=test_df,
    target_col=target_col,
    task_type="regression",
    timestamp=timestamp,
    # folds=(ensemble_size, slurm_array_task_id),
    epochs=epochs,
    checkpoint="wandb",  # None | 'local' | 'wandb',
    input_col=input_col,
    learning_rate=learning_rate,
    batch_size=batch_size,
    wandb_path=WANDB_PATH,
    run_params=run_params,
)
