"""Train a Wrenformer ensemble on target_col of data_path."""

# %%
import os
from importlib.metadata import version

import pandas as pd
from aviary.train import df_train_test_split, train_wrenformer
from pymatviz.enums import Key

from matbench_discovery import WANDB_PATH, timestamp, today
from matbench_discovery.enums import DataFiles, MbdKey, Model
from matbench_discovery.hpc import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-08-13"


# %%
epochs = 300
data_path = DataFiles.mp_energies.path
target_col = Key.form_energy
# data_path = f"{ROOT}/data/2022-08-25-m3gnet-trainset-mp-2021-struct-energy.json.gz"
# target_col = "mp_energy_per_atom"
data_name = "m3gnet-trainset" if "m3gnet" in data_path else "mp"
ensemble_size = 10
job_name = f"{today}-train-{Model.wrenformer}-ens={ensemble_size}-robust-{data_name}"
dataset = "mp"
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
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))

print(f"\nJob {job_name} started {timestamp}")

df_data = pd.read_csv(data_path).set_index(Key.mat_id, drop=False)

if target_col not in df_data:
    raise TypeError(f"{target_col!s} not in {df_data.columns=}")
if MbdKey.wyckoff_spglib not in df_data:
    raise TypeError(f"{MbdKey.wyckoff_spglib!s} not in {df_data.columns=}")
train_df, test_df = df_train_test_split(df_data, test_size=0.05)

run_params = dict(
    data_path=data_path,
    versions={dep: version(dep) for dep in ("aviary", "numpy", "torch")},
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
    input_col=MbdKey.wyckoff_spglib,
    learning_rate=learning_rate,
    batch_size=batch_size,
    wandb_path=WANDB_PATH,
    run_params=run_params,
)
