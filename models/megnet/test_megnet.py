"""Get MEGNet formation energy predictions on WBM test set.
To slurm submit this file: python path/to/file.py slurm-submit
Requires MEGNet installation: pip install megnet
See https://github.com/materialsvirtuallab/megnet.
"""

# %%
import os
from importlib.metadata import version

import numpy as np
import pandas as pd
import pymatviz as pmv
import wandb
from megnet.utils.models import load_model
from pymatgen.core import Structure
from pymatviz.enums import Key
from sklearn.metrics import r2_score
from tqdm import tqdm

from matbench_discovery import timestamp, today
from matbench_discovery.data import DataFiles, Model, df_wbm
from matbench_discovery.enums import MbdKey, Task
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-11-14"

task_type = Task.RS2RE
module_dir = os.path.dirname(__file__)
job_name = f"megnet-wbm-{task_type}"
out_path = os.getenv("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}.csv.gz")
slurm_array_task_count = 1

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exciting early")

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=module_dir,
    account="matgen",
    time="11:55:0",
    slurm_flags="--mem 30G",
    array=f"1-{slurm_array_task_count}",
    # TF_CPP_MIN_LOG_LEVEL=2 means INFO and WARNING logs are not printed
    # https://stackoverflow.com/a/40982782
    pre_cmd="TF_CPP_MIN_LOG_LEVEL=2",
)


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
data_path = {
    Task.IS2RE: DataFiles.wbm_initial_structures.path,
    Task.RS2RE: DataFiles.wbm_computed_structure_entries.path,
    "chgnet_structure": Model.chgnet.path.replace(".csv.gz", ".json.gz"),
    "m3gnet_structure": Model.m3gnet.path.replace(".csv.gz", ".json.gz"),
}[task_type]
print(f"\nJob {job_name} started {timestamp}")
print(f"{data_path=}")
if MbdKey.e_form_dft not in df_wbm:
    raise KeyError(f"{MbdKey.e_form_dft!s} not in {df_wbm.columns=}")

df_in = pd.read_json(data_path).set_index(Key.mat_id)
if slurm_array_task_count > 1:
    df_in = np.array_split(df_in, slurm_array_task_count)[slurm_array_task_id - 1]
megnet_mp_e_form = load_model(model_name := "Eform_MP_2019")


# %%
run_params = {
    "data_path": data_path,
    "versions": {dep: version(dep) for dep in ("megnet", "numpy")},
    "model_name": model_name,
    Key.task_type: task_type,
    "target_col": MbdKey.e_form_dft,
    "df": {"shape": str(df_in.shape), "columns": ", ".join(df_in)},
    "slurm_vars": slurm_vars,
    Key.model_params: sum(
        np.prod(p.shape) for p in megnet_mp_e_form.trainable_variables
    ),
}

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %% input_col=task_type for CHGNet and M3GNet
input_col = {Task.IS2RE: Key.init_struct, Task.RS2RE: Key.final_struct}.get(
    task_type, task_type
)

if task_type == Task.RS2RE:
    df_in[input_col] = [cse["structure"] for cse in df_in[Key.cse]]

structures = df_in[input_col].map(Structure.from_dict).to_dict()

megnet_e_form_preds = {}
for material_id in tqdm(structures):
    if material_id in megnet_e_form_preds:
        continue
    try:
        structure = structures[material_id]
        e_form_per_atom = megnet_mp_e_form.predict_structure(structure)[0]
        megnet_e_form_preds[material_id] = e_form_per_atom
    except Exception as exc:
        print(f"Failed to predict {material_id=}: {exc}")


# %%
print(f"{len(megnet_e_form_preds)=:,}")
print(f"{len(structures)=:,}")
print(f"missing: {len(structures) - len(megnet_e_form_preds):,}")
pred_col = "e_form_per_atom_megnet"
# remove legacy MP corrections that MEGNet was trained on and apply newer MP2020
# corrections instead
df_megnet = pd.DataFrame({f"{pred_col}_old_corr": megnet_e_form_preds})
df_megnet[pred_col] = (
    pd.Series(megnet_e_form_preds)
    - df_wbm.e_correction_per_atom_mp_legacy
    + df_wbm.e_correction_per_atom_mp2020
)
df_megnet.index.name = Key.mat_id
if task_type != Task.IS2RE:
    df_megnet = df_megnet.add_suffix(f"_{task_type.lower()}")


df_megnet.add_suffix(f"_{task_type.lower()}").round(4).to_csv(out_path)

# df_megnet = pd.read_csv(
#     f"{ROOT}/models/{Model.megnet.path}"
# ).set_index(Key.mat_id)


# %% compare MEGNet predictions with old and new MP corrections
ax = pmv.density_scatter(df=df_megnet, x=pred_col, y=f"{pred_col}_old_corr")
pmv.save_fig(ax, "megnet-e-form-preds-old-vs-new-corr.png")


# %%
df_wbm[pred_col] = df_megnet[pred_col]
table = wandb.Table(dataframe=df_wbm[[MbdKey.e_form_dft, pred_col]].reset_index())

MAE = (df_wbm[MbdKey.e_form_dft] - df_wbm[pred_col]).abs().mean()
R2 = r2_score(*df_wbm[[MbdKey.e_form_dft, pred_col]].dropna().to_numpy().T)
title = f"{model_name} {task_type} {MAE=:.4} {R2=:.4}"
print(title)

wandb_scatter(table, fields=dict(x=MbdKey.e_form_dft, y=pred_col), title=title)
