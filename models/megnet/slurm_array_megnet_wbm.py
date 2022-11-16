# %%
from __future__ import annotations

import os
from datetime import datetime
from importlib.metadata import version

import pandas as pd
import wandb
from megnet.utils.models import load_model
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.plot_scripts import df_wbm
from matbench_discovery.slurm import slurm_submit

"""
To slurm submit this file: python path/to/file.py slurm-submit
Requires Megnet installation: pip install megnet
https://github.com/materialsvirtuallab/megnet
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-11-14"

task_type = "IS2RE"
timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@")[0]
module_dir = os.path.dirname(__file__)
job_name = f"megnet-wbm-{task_type}"
out_dir = f"{module_dir}/{today}-{job_name}"

slurm_vars = slurm_submit(
    job_name=job_name,
    log_dir=out_dir,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time=(slurm_max_job_time := "12:0:0"),
    # TF_CPP_MIN_LOG_LEVEL=2 means INFO and WARNING logs are not printed
    # https://stackoverflow.com/a/40982782
    pre_cmd="TF_CPP_MIN_LOG_LEVEL=2",
)


# %%
print(f"Job started running {timestamp}")

out_path = f"{out_dir}/megnet-e-form-preds.csv"
if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")


# %%
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
print(f"Loading from {data_path=}")
df_wbm_structs = pd.read_json(data_path).set_index("material_id")


megnet_mp_e_form = load_model(model_name := "Eform_MP_2019")


# %%
run_params = dict(
    data_path=data_path,
    megnet_version=version("megnet"),
    model_name=model_name,
    task_type=task_type,
    slurm_max_job_time=slurm_max_job_time,
    df=dict(shape=str(df_wbm_structs.shape), columns=", ".join(df_wbm_structs)),
    slurm_vars=slurm_vars,
)
if wandb.run is None:
    wandb.login()

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
for material_id, structure in tqdm(structures.items(), total=len(structures)):
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
out_col = "e_form_per_atom_megnet"
df_wbm[out_col] = pd.Series(megnet_e_form_preds)

df_wbm[out_col].reset_index().to_csv(out_path)


# %%
fields = {"x": "e_form_per_atom_mp2020_corrected", "y": out_col}
cols = list(fields.values())
assert all(col in df_wbm for col in cols)

table = wandb.Table(dataframe=df_wbm[cols].reset_index())

MAE = (df_wbm[fields["x"]] - df_wbm[fields["y"]]).abs().mean()

scatter_plot = wandb.plot_table(
    vega_spec_name="janosh/scatter-parity",
    data_table=table,
    fields=fields,
    string_fields={"title": f"{model_name} {task_type} {MAE=:.4}"},
)

wandb.log({"true_pred_scatter": scatter_plot})
