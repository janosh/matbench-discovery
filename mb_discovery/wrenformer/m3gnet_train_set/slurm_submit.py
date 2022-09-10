# %%
import os
import subprocess
from datetime import datetime

from aviary import ROOT

"""
Write a Python job file and sbatch it using subprocess.run()
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-08-25"


# %%
epochs = 300
n_attn_layers = 3
embedding_aggregations = ("mean",)
optimizer = "AdamW"
lr = 3e-4
n_folds = 10
df_or_path = f"{ROOT}/datasets/2022-08-25-m3gnet-trainset-mp-2021-struct-energy.json.gz"
target_col = "mp_energy_per_atom"
task_type = "regression"
checkpoint = "wandb"  # None | 'local' | 'wandb'
batch_size = 128
model_name = f"wrenformer-robust-{epochs=}"
swa_start = None


log_dir = f"{os.path.dirname(__file__)}/job-logs"
os.makedirs(log_dir, exist_ok=True)
timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"

python_script = f"""import os
from datetime import datetime

from aviary.wrenformer.train import train_wrenformer_on_df

print(f"Job started running {{datetime.now():%Y-%m-%d@%H-%M}}")
job_id = os.environ["SLURM_JOB_ID"]
print(f"{{job_id=}}")
print("{model_name=}")
print("{df_or_path=}")

job_array_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
print(f"{{job_array_id=}}")

train_wrenformer_on_df(
    {model_name=},
    {target_col=},
    {df_or_path=},
    {timestamp=},
    test_size=0.05,
    # folds=(n_folds, job_array_id),
    {epochs=},
    {n_attn_layers=},
    {checkpoint=},
    {optimizer=},
    learning_rate={lr},
    {embedding_aggregations=},
    {batch_size=},
    {swa_start=},
    wandb_path="aviary/wrenformer-on-m3gnet-trainset",
)
"""


submit_script = f"{log_dir}/{timestamp}-{model_name}.py"

# prepend into sbatch script to source module command and load default env
# for Ampere GPU partition before actual job command
slurm_setup = ". /etc/profile.d/modules.sh; module load rhel8/default-amp;"


# %%
slurm_cmd = f"""sbatch
    --partition ampere
    --account LEE-SL3-GPU
    --time 8:0:0
    --nodes 1
    --gpus-per-node 1
    --chdir {log_dir}
    --array 0-{n_folds - 1}
    --out {timestamp}-{model_name}-%A-%a.log
    --job-name {model_name}
    --wrap '{slurm_setup} python {submit_script}'
"""

header = f'"""\nSlurm submission command:\n{slurm_cmd}"""'

with open(submit_script, "w") as file:
    file.write(f"{header}\n\n{python_script}")


# %% submit slurm job
result = subprocess.run(slurm_cmd.replace("\n    ", " "), check=True, shell=True)
print(f"\n{submit_script = }")
