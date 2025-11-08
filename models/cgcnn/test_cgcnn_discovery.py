"""
Download WandB checkpoints for an ensemble of CGCNN models trained on all MP
formation energies, then make predictions on some dataset, prints ensemble metrics and
saves predictions to CSV.
"""

# %%
import os
from importlib.metadata import version

import pandas as pd
import torch
import wandb
from aviary.cgcnn.data import CrystalGraphData, collate_batch
from aviary.cgcnn.model import CrystalGraphConvNet
from aviary.predict import predict_from_wandb_checkpoints
from pymatgen.core import Structure
from pymatviz.enums import Key
from torch.utils.data import DataLoader
from tqdm import tqdm

from matbench_discovery import CHECKPOINT_DIR, WANDB_PATH, WBM_DIR, today
from matbench_discovery.data import df_wbm
from matbench_discovery.enums import DataFiles, MbdKey, Model, Task
from matbench_discovery.hpc import slurm_submit
from matbench_discovery.plots import wandb_scatter

__author__ = "Janosh Riebesell"
__date__ = "2022-08-15"


task_type = Task.IS2RE
debug = False
model_name = Model.cgcnn  # or Model.cgcnn_p
job_name = f"{model_name}/{today}-wbm-{task_type}"
module_dir = os.path.dirname(__file__)
out_dir = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    account="matgen",
    time="2:0:0",
    out_dir=out_dir,
    slurm_flags="--nodes 1 --gpus-per-node 1",
)


# %%
data_path = {
    Task.IS2RE: DataFiles.wbm_initial_structures.path,
    Task.RS2RE: DataFiles.wbm_computed_structure_entries.path,
    "IS2RE-debug": f"{WBM_DIR}/2022-10-19-wbm-init-structs.json-1k-samples.bz2",
}[task_type]
input_col = {Task.IS2RE: Key.initial_struct, Task.RS2RE: Key.final_struct}[task_type]

df_in = pd.read_json(data_path, lines=True).set_index(Key.mat_id)
if input_col not in df_in:
    raise TypeError(f"{input_col!s} not in {df_in.columns=}")

df_in[MbdKey.e_form_dft] = df_wbm[MbdKey.e_form_dft]
if task_type == Task.RS2RE:
    df_in[input_col] = [cse["structure"] for cse in df_in[Key.computed_structure_entry]]

df_in[input_col] = [
    Structure.from_dict(dct) for dct in tqdm(df_in[input_col], disable=None)
]

filters = {
    # "display_name": {"$regex": "^train-cgcnn-augment=3-"},
    # "created_at": {"$gt": "2022-12-03", "$lt": "2022-12-04"},
    "display_name": {"$regex": "^train-cgcnn-augment=0-"},
    "created_at": {"$gt": "2023-01-09", "$lt": "2023-01-10"},
}
runs = wandb.Api().runs(WANDB_PATH, filters=filters)
expected_runs = 10
if len(runs) != expected_runs:
    raise ValueError(
        f"{expected_runs=}, got {len(runs)} filtering {WANDB_PATH=} with {filters=}"
    )

for idx, run in enumerate(runs):
    for key, val in run.config.items():
        if val == runs[0].config[key] or key.startswith(("slurm_", "timestamp")):
            continue
        raise ValueError(
            f"Run configs not identical: runs[{idx}][{key}]={val}, {runs[0][key]=}"
        )

run_params = dict(
    data_path=data_path,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    versions={dep: version(dep) for dep in ("aviary", "numpy", "torch")},
    ensemble_size=len(runs),
    task_type=task_type,
    target_col=MbdKey.e_form_dft,
    input_col=input_col,
    wandb_run_filters=filters,
    slurm_vars=slurm_vars,
    training_run_ids=[run.id for run in runs],
)
try:  # load checkpoint to get number of parameters
    runs[0].file("checkpoint.pth").download(root=module_dir)
    state_dict = torch.load(f"{module_dir}/checkpoint.pth", map_location="cpu")
    model = CrystalGraphConvNet(**state_dict["model_params"])
    run_params[Key.model_params] = model.num_params
except Exception as exc:
    print(exc)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)

cg_data = CrystalGraphData(
    df_in, task_dict={MbdKey.e_form_dft: "regression"}, structure_col=input_col
)
data_loader = DataLoader(
    cg_data, batch_size=1024, shuffle=False, collate_fn=collate_batch
)


# %%
out = predict_from_wandb_checkpoints(
    runs,
    # dropping isolated-atom structs means len(cg_data.df) < len(df)
    cache_dir=CHECKPOINT_DIR,
    df=cg_data.df.drop(columns=input_col),
    target_col=MbdKey.e_form_dft,
    model_cls=CrystalGraphConvNet,
    data_loader=data_loader,
)
# type narrow for ty's benefit, better would be to properly overload predict_from_wandb_checkpoints  # noqa: E501
if not (isinstance(out, tuple) and len(out) == 2):
    raise TypeError(f"{out=} should be 2-tuple")
df_in, ensemble_metrics = out

slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")
df_in.round(4).to_csv(f"{out_dir}/{job_name}-preds-{slurm_array_job_id}.csv.gz")
pred_col = f"{MbdKey.e_form_dft}_pred_ens"
if pred_col not in df_in:
    raise KeyError(f"{pred_col} not in {df_in.columns=}")
table = wandb.Table(dataframe=df_in[[MbdKey.e_form_dft, pred_col]].reset_index())


# %%
MAE = ensemble_metrics.MAE.mean()
R2 = ensemble_metrics.R2.mean()

title = f"CGCNN {task_type} ensemble={len(runs)} {MAE=:.4} {R2=:.4}"

wandb_scatter(table, fields=dict(x=MbdKey.e_form_dft, y=pred_col), title=title)
