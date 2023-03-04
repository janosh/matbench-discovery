"""Get MEGNet formation energy predictions on WBM test set.
To slurm submit this file: python path/to/file.py slurm-submit
Requires MEGNet installation: pip install megnet
See https://github.com/materialsvirtuallab/megnet.
"""


# %%
from __future__ import annotations

import os
from importlib.metadata import version

import pandas as pd
import wandb
from megnet.utils.models import load_model
from pymatgen.core import Structure
from sklearn.metrics import r2_score
from tqdm import tqdm

from matbench_discovery import DEBUG, timestamp, today
from matbench_discovery.data import DATA_FILES, df_wbm
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2022-11-14"

task_type = "IS2RE"
module_dir = os.path.dirname(__file__)
job_name = f"megnet-wbm-{task_type}{'-debug' if DEBUG else ''}"
out_dir = os.environ.get("SBATCH_OUTPUT", f"{module_dir}/{today}-{job_name}")

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time="12:0:0",
    slurm_flags=("--mem", "30G"),
    # TF_CPP_MIN_LOG_LEVEL=2 means INFO and WARNING logs are not printed
    # https://stackoverflow.com/a/40982782
    pre_cmd="TF_CPP_MIN_LOG_LEVEL=2",
)


# %%
out_path = f"{out_dir}/megnet-e-form-preds.csv"
if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")

data_path = {
    "IS2RE": DATA_FILES.wbm_initial_structures,
    "RS2RE": DATA_FILES.wbm_computed_structure_entries,
}[task_type]
print(f"\nJob started running {timestamp}")
print(f"{data_path=}")
e_form_col = "e_form_per_atom_mp2020_corrected"
assert e_form_col in df_wbm, f"{e_form_col=} not in {list(df_wbm)=}"

df_in = pd.read_json(data_path).set_index("material_id")
megnet_mp_e_form = load_model(model_name := "Eform_MP_2019")


# %%
run_params = dict(
    data_path=data_path,
    megnet_version=version("megnet"),
    numpy_version=version("numpy"),
    model_name=model_name,
    task_type=task_type,
    target_col=e_form_col,
    df=dict(shape=str(df_in.shape), columns=", ".join(df_in)),
    slurm_vars=slurm_vars,
)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %%
input_col = {"IS2RE": "initial_structure", "RS2RE": "relaxed_structure"}[task_type]

if task_type == "RS2RE":
    df_in[input_col] = [x["structure"] for x in df_in.computed_structure_entry]

structures = df_in[input_col].map(Structure.from_dict).to_dict()

megnet_e_form_preds = {}
for material_id in tqdm(structures, disable=None):
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
# old columns contains direct MEGNet predictions which was trained on legacy-corrected
# MP formation energies
df_wbm[f"{pred_col}_old"] = pd.Series(megnet_e_form_preds)

# remove legacy MP corrections that MEGNet was trained on and apply newer MP2020
# corrections instead
df_wbm[pred_col] = (
    df_wbm[pred_col]
    - df_wbm.e_correction_per_atom_mp_legacy
    + df_wbm.e_correction_per_atom_mp2020
)

df_wbm.filter(like=pred_col).round(4).to_csv(
    "2022-11-18-megnet-wbm-IS2RE/megnet-e-form-preds.csv"
)

# df_megnet = pd.read_csv(f"{ROOT}/models/{PRED_FILES.megnet}").set_index("material_id")


# %%
table = wandb.Table(dataframe=df_wbm[[e_form_col, pred_col]].reset_index())

MAE = (df_wbm[e_form_col] - df_wbm[pred_col]).abs().mean()
R2 = r2_score(df_wbm[e_form_col], df_wbm[pred_col])
title = f"{model_name} {task_type} {MAE=:.4} {R2=:.4}"
print(title)

wandb_scatter(table, fields=dict(x=e_form_col, y=pred_col), title=title)
