# %%
import os
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

today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)


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
data_path = f"{ROOT}/data/2022-09-16-mp-computed-structure-entries.json.gz"
# data_path = f"{ROOT}/data/wbm/2022-10-19-wbm-cses+init-structs.json.bz2"
df = pd.read_json(data_path).set_index("material_id")

df["structure"] = [Structure.from_dict(x["structure"]) for x in tqdm(df.entry)]


# %%
df_featurized = featurizer.featurize_dataframe(df, "structure", ignore_errors=True)


# %%
df_featurized.to_json(f"{module_dir}/{today}mp-train-voronoi-tesselation.json.gz")
