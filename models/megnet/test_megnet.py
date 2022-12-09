# %%
from __future__ import annotations

import os
from importlib.metadata import version

import pandas as pd
import wandb
from megnet.utils.models import load_model
from sklearn.metrics import r2_score
from tqdm import tqdm

from matbench_discovery import DEBUG, ROOT, timestamp, today
from matbench_discovery.data import df_wbm
from matbench_discovery.plots import wandb_scatter
from matbench_discovery.slurm import slurm_submit

"""
To slurm submit this file: python path/to/file.py slurm-submit
Requires MEGNet installation: pip install megnet
https://github.com/materialsvirtuallab/megnet
"""

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

data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
print(f"\nJob started running {timestamp}")
print(f"{data_path=}")
target_col = "e_form_per_atom_mp2020_corrected"
assert target_col in df_wbm, f"{target_col=} not in {list(df_wbm)=}"

df_wbm_structs = pd.read_json(data_path).set_index("material_id")
megnet_mp_e_form = load_model(model_name := "Eform_MP_2019")


# %%
run_params = dict(
    data_path=data_path,
    megnet_version=version("megnet"),
    model_name=model_name,
    task_type=task_type,
    target_col=target_col,
    df=dict(shape=str(df_wbm_structs.shape), columns=", ".join(df_wbm_structs)),
    slurm_vars=slurm_vars,
)

wandb.init(project="matbench-discovery", name=job_name, config=run_params)


# %%
if task_type == "IS2RE":
    from pymatgen.core import Structure

    structures = df_wbm_structs.initial_structure.map(Structure.from_dict)
elif task_type == "RS2RE":
    from pymatgen.entries.computed_entries import ComputedStructureEntry

    df_wbm_structs.cse = df_wbm_structs.cse.map(ComputedStructureEntry.from_dict)
    structures = df_wbm_structs.cse.map(lambda x: x.structure)
else:
    raise ValueError(f"Unknown {task_type = }")

megnet_e_form_preds = {}
for material_id, structure in tqdm(
    structures.items(), disable=None, total=len(structures)
):
    if material_id in megnet_e_form_preds:
        continue
    try:
        e_form_per_atom = megnet_mp_e_form.predict_structure(structure)[0]
        megnet_e_form_preds[material_id] = e_form_per_atom
    except Exception as exc:
        print(f"Failed to predict {material_id=}: {exc}")


# %%
print(f"{len(megnet_e_form_preds)=:,}")
print(f"{len(structures)=:,}")
print(f"missing: {len(structures) - len(megnet_e_form_preds):,}")
pred_col = "e_form_per_atom_megnet"
df_wbm[pred_col] = pd.Series(megnet_e_form_preds)

df_wbm[pred_col].reset_index().round(4).to_csv(out_path, index=False)


# %%
table = wandb.Table(dataframe=df_wbm[[target_col, pred_col]].reset_index())

MAE = (df_wbm[target_col] - df_wbm[pred_col]).abs().mean()
R2 = r2_score(df_wbm[target_col], df_wbm[pred_col])
title = f"{model_name} {task_type} {MAE=:.4} {R2=:.4}"
print(title)

wandb_scatter(table, fields=dict(x=target_col, y=pred_col), title=title)
