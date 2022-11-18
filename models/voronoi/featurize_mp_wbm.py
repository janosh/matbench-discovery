# %%
import os
import warnings
from datetime import datetime

import matminer.featurizers.composition as feat_comp
import matminer.featurizers.structure as feat_struct
import numpy as np
import pandas as pd
import wandb
from matminer.featurizers.base import MultipleFeaturizer
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import ROOT, as_dict_handler
from matbench_discovery.slurm import slurm_submit

today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)


# data_path = f"{ROOT}/data/mp/2022-09-16-mp-computed-structure-entries.json.gz"
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
input_col = "initial_structure"
data_name = "wbm" if "wbm" in data_path else "mp"
slurm_array_task_count = 10
job_name = f"voronoi-features-{data_name}"
log_dir = f"{module_dir}/{today}-{job_name}"

slurm_vars = slurm_submit(
    job_name=job_name,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time=(slurm_max_job_time := "5:0:0"),
    array=f"1-{slurm_array_task_count}",
    log_dir=log_dir,
)


# %%
slurm_array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
run_name = f"{job_name}-{slurm_array_task_id}"
out_path = f"{log_dir}/{run_name}.csv.bz2"

if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")

df = pd.read_json(data_path).set_index("material_id")
df_this_job: pd.DataFrame = np.array_split(df, slurm_array_task_count)[
    slurm_array_task_id - 1
]

if data_name == "mp":
    struct_dicts = [x["structure"] for x in df_this_job.entry]
if data_name == "wbm":
    struct_dicts = df_this_job.initial_structure

df_this_job[input_col] = [
    Structure.from_dict(x) for x in tqdm(df_this_job.initial_structure, disable=None)
]


# %%
run_params = dict(
    data_path=data_path,
    slurm_max_job_time=slurm_max_job_time,
    df=dict(shape=str(df_this_job.shape), columns=", ".join(df_this_job)),
    input_col=input_col,
    slurm_vars=slurm_vars,
)
if wandb.run is None:
    wandb.login()

wandb.init(
    project="matbench-discovery",
    name=run_name,
    config=run_params,
)


# %% Create the featurizer: Ward et al. use a variety of different featurizers
# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.024104
featurizers = [
    feat_struct.SiteStatsFingerprint.from_preset("CoordinationNumber_ward-prb-2017"),
    feat_struct.StructuralHeterogeneity(),
    feat_struct.ChemicalOrdering(),
    feat_struct.MaximumPackingEfficiency(),
    feat_struct.SiteStatsFingerprint.from_preset(
        "LocalPropertyDifference_ward-prb-2017"
    ),
    feat_struct.StructureComposition(feat_comp.Stoichiometry()),
    feat_struct.StructureComposition(feat_comp.ElementProperty.from_preset("magpie")),
    feat_struct.StructureComposition(feat_comp.ValenceOrbital(props=["frac"])),
    feat_struct.StructureComposition(feat_comp.IonProperty(fast=True)),
]
featurizer = MultipleFeaturizer(featurizers)
# multiprocessing seems to be the cause of OOM errors on large structures even when
# taking only small slice of the data and launching slurm jobs with --mem 100G
featurizer.set_n_jobs(1)


# %% prints lots of pymatgen warnings
# > No electronegativity for Ne. Setting to NaN. This has no physical meaning, ...
warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")

df_features = featurizer.featurize_dataframe(
    df_this_job, input_col, ignore_errors=True, pbar=dict(position=0, leave=True)
).drop(columns=input_col)


# %%
df_features.to_csv(out_path, default_handler=as_dict_handler)

wandb.log({"voronoi_features": wandb.Table(dataframe=df_features)})
