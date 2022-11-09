# %%
import os
import warnings
from datetime import datetime

import pandas as pd
from matminer.featurizers.base import MultipleFeaturizer
from matminer.featurizers.composition import (
    ElementProperty,
    IonProperty,
    Stoichiometry,
    ValenceOrbital,
)
from matminer.featurizers.structure import (
    ChemicalOrdering,
    MaximumPackingEfficiency,
    SiteStatsFingerprint,
    StructuralHeterogeneity,
    StructureComposition,
)
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.slurm import slurm_submit_python

today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)


# data_path = f"{ROOT}/data/mp/2022-09-16-mp-computed-structure-entries.json.gz"
data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-init-structs.json.bz2"
input_col = "structure"
data_name = "wbm" if "wbm" in data_path else "mp"

slurm_submit_python(
    job_name=f"voronoi-featurize-{data_name}",
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time="3:0:0",
    log_dir=module_dir,
    slurm_flags=("--mem=40G",),
)


# %% Create the featurizer: Ward et al. use a variety of different featurizers
# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.024104
featurizer = MultipleFeaturizer(
    [
        SiteStatsFingerprint.from_preset("CoordinationNumber_ward-prb-2017"),
        StructuralHeterogeneity(),
        ChemicalOrdering(),
        MaximumPackingEfficiency(),
        SiteStatsFingerprint.from_preset("LocalPropertyDifference_ward-prb-2017"),
        StructureComposition(Stoichiometry()),
        StructureComposition(ElementProperty.from_preset("magpie")),
        StructureComposition(ValenceOrbital(props=["frac"])),
        StructureComposition(IonProperty(fast=True)),
    ],
)


# %%
df = pd.read_json(data_path).set_index("material_id")

if data_name == "mp":
    struct_dicts = [x["structure"] for x in df.entry]
if data_name == "wbm":
    struct_dicts = df.initial_structure

df[input_col] = [
    Structure.from_dict(x) for x in tqdm(df.initial_structure, disable=None)
]


# %% prints lots of pymatgen warnings
# > No electronegativity for Ne. Setting to NaN. This has no physical meaning, ...
warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")

df_features = featurizer.featurize_dataframe(
    df, input_col, ignore_errors=True, pbar=True
)


# %%
df_features.to_json(
    f"{module_dir}/{today}-voronoi-tesselation-{data_name}-features.json.gz"
)
