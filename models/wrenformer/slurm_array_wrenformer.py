# %%
import os
from datetime import datetime

from aviary.wrenformer.train import train_wrenformer_on_df

from matbench_discovery import ROOT
from matbench_discovery.slurm import slurm_submit_python

"""
Train a Wrenformer
ensemble of size n_folds on target_col of df_or_path.
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-08-13"


# %%
df_or_path = f"{ROOT}/data/2022-08-13-mp-energies.json.gz"
target_col = "energy_per_atom"
# df_or_path = f"{ROOT}/data/2022-08-25-m3gnet-trainset-mp-2021-struct-energy.json.gz"
# target_col = "mp_energy_per_atom"

epochs = 300
job_name = f"wrenformer-robust-{epochs=}-{target_col}"
n_folds = 10
today = f"{datetime.now():%Y-%m-%d}"
dataset = "mp"
# dataset = 'm3gnet_train_set'
log_dir = f"{os.path.dirname(__file__)}/{dataset}/{today}-{job_name}"

slurm_submit_python(
    job_name=job_name,
    partition="ampere",
    time="8:0:0",
    array=f"1-{n_folds}",
    log_dir=log_dir,
    account="LEE-SL3-GPU",
    slurm_flags=("--nodes 1", "--gpus-per-node 1"),
    # prepend into sbatch script to source module command and load default env
    # for Ampere GPU partition before actual job command
    pre_cmd=". /etc/profile.d/modules.sh; module load rhel8/default-amp;",
    # if running on CPU, unsetting OMP threads allows using PyTorch to use all cores
    # https://docs.hpc.cam.ac.uk/hpc/software-packages/pytorch.html
    # pre_cmd="unset OMP_NUM_THREADS",
)


# %%
n_attn_layers = 3
embedding_aggregations = ("mean",)
optimizer = "AdamW"
learning_rate = 3e-4
task_type = "regression"
checkpoint = "wandb"  # None | 'local' | 'wandb'
batch_size = 128
swa_start = None
timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"

print(f"Job started running {datetime.now():%Y-%m-%d@%H-%M}")
slurm_job_id = os.environ.get("SLURM_JOB_ID")
slurm_array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))

print(f"{slurm_job_id=}")
print(f"{slurm_array_task_id=}")
print(f"{job_name=}")
print(f"{df_or_path=}")

train_wrenformer_on_df(
    run_name=job_name,
    target_col=target_col,
    df_or_path=df_or_path,
    timestamp=timestamp,
    test_size=0.05,
    # folds=(n_folds, slurm_array_task_id),
    epochs=epochs,
    n_attn_layers=n_attn_layers,
    checkpoint=checkpoint,
    optimizer=optimizer,
    learning_rate=learning_rate,
    embedding_aggregations=embedding_aggregations,
    batch_size=batch_size,
    swa_start=swa_start,
    wandb_path="janosh/matbench-discovery",
)
