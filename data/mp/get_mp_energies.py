# %%
import os

import pandas as pd
from aviary.utils import as_dict_handler
from aviary.wren.utils import get_aflow_label_from_spglib
from mp_api.client import MPRester
from tqdm import tqdm

from matbench_discovery import today

"""
Download all MP formation and above hull energies on 2022-08-13.

Related EDA of MP formation energies:
https://github.com/janosh/pymatviz/blob/main/examples/mp_bimodal_e_form.ipynb
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-08-13"


module_dir = os.path.dirname(__file__)


# %% query all MP formation energies on 2022-08-13
fields = [
    "material_id",
    "task_ids",
    "formula_pretty",
    "formation_energy_per_atom",
    "energy_per_atom",
    "structure",
    "symmetry",
    "energy_above_hull",
]
with MPRester(use_document_model=False) as mpr:
    docs = mpr.summary.search(fields=fields)

print(f"{today}: {len(docs) = :,}")
# 2022-08-13: len(docs) = 146,323


# %%
df = pd.DataFrame(docs).set_index("material_id")
df.pop("_id")

df["spacegroup_number"] = df.pop("symmetry").map(lambda x: x.number)

df["wyckoff_spglib"] = [get_aflow_label_from_spglib(x) for x in tqdm(df.structure)]

df.to_json(f"{module_dir}/{today}-mp-energies.json.gz", default_handler=as_dict_handler)

# df = pd.read_json(f"{module_dir}/2022-08-13-mp-energies.json.gz")