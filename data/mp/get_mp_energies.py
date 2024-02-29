# %%
import os

import pandas as pd
from aviary.wren.utils import get_aflow_label_from_spglib
from mp_api.client import MPRester
from pymatgen.core import Structure
from pymatviz.utils import annotate_metrics
from tqdm import tqdm

from matbench_discovery import STABILITY_THRESHOLD, today
from matbench_discovery.data import DATA_FILES
from matbench_discovery.enums import Key

"""
Download all MP formation and above hull energies on 2023-01-10.

Related EDA of MP formation energies:
https://github.com/janosh/pymatviz/blob/-/examples/mp_bimodal_e_form.ipynb
"""

__author__ = "Janosh Riebesell"
__date__ = "2023-01-10"

module_dir = os.path.dirname(__file__)


# %%
fields = {
    Key.mat_id,
    "formula_pretty",
    Key.form_energy,
    "energy_per_atom",
    "symmetry",
    "energy_above_hull",
    "decomposition_enthalpy",
    "energy_type",
    "nsites",
}

with MPRester(use_document_model=False) as mpr:
    docs = mpr.thermo.search(fields=fields, thermo_types=["GGA_GGA+U"])

assert fields == set(docs[0]), f"missing fields: {fields - set(docs[0])}"
print(f"{today}: {len(docs) = :,}")
# 2022-08-13: len(docs) = 146,323
# 2023-01-10: len(docs) = 154,718


# %%
df = pd.DataFrame(docs).set_index(Key.mat_id)
df = df.rename(columns={"formula_pretty": Key.formula, "nsites": Key.n_sites})

df_spg = pd.json_normalize(df.pop("symmetry"))[["number", "symbol"]]
df["spacegroup_symbol"] = df_spg.symbol.to_numpy()

df.energy_type.value_counts().plot.pie(backend="plotly", autopct="%1.1f%%")
# GGA: 72.2%, GGA+U: 27.8%


# %%
df_cse = pd.read_json(DATA_FILES.mp_computed_structure_entries).set_index(Key.mat_id)

df_cse[Key.struct] = [
    Structure.from_dict(cse[Key.struct]) for cse in tqdm(df_cse.entry)
]
df_cse[Key.wyckoff] = [
    get_aflow_label_from_spglib(struct, errors="ignore")
    for struct in tqdm(df_cse.structure)
]
# make sure symmetry detection succeeded for all structures
assert df_cse[Key.wyckoff].str.startswith("invalid").sum() == 0
df[Key.wyckoff] = df_cse[Key.wyckoff]

spg_nums = df[Key.wyckoff].str.split("_").str[2].astype(int)
# make sure all our spacegroup numbers match MP's
assert (spg_nums.sort_index() == df_spg["number"].sort_index()).all()

df.to_csv(DATA_FILES.mp_energies)
# df = pd.read_csv(DATA_FILES.mp_energies, na_filter=False).set_index(Key.mat_id)


# %% reproduce fig. 1b from https://arxiv.org/abs/2001.10591 (as data consistency check)
ax = df.plot.scatter(
    x=Key.form_energy,
    y="decomposition_enthalpy",
    alpha=0.1,
    xlim=[-5, 1],
    ylim=[-1, 1],
    color=(df.decomposition_enthalpy > STABILITY_THRESHOLD).map(
        {True: "red", False: "blue"}
    ),
    title=f"{today} - {len(df):,} MP entries",
)

annotate_metrics(df.formation_energy_per_atom, df.decomposition_enthalpy)
# result on 2023-01-10: plots match. no correlation between formation energy and
# decomposition enthalpy. R^2 = -1.571, MAE = 1.604
# ax.figure.savefig(f"{module_dir}/mp-decomp-enth-vs-e-form.webp", dpi=300)


# %% scatter plot energy above convex hull vs decomposition enthalpy
# https://berkeleytheory.slack.com/archives/C16RE1TUN/p1673887564955539
mask_above_line = df.energy_above_hull - df.decomposition_enthalpy.clip(0) > 0.1
ax = df.plot.scatter(
    x="decomposition_enthalpy",
    y="energy_above_hull",
    color=mask_above_line.map({True: "red", False: "blue"}),
    hover_data=["index", Key.formula, Key.form_energy],
)
# most points lie on line y=x for x > 0 and y = 0 for x < 0.
n_above_line = sum(mask_above_line)
ax.set(
    title=f"{n_above_line:,} / {len(df):,} = {n_above_line / len(df):.1%} "
    "MP materials with\nenergy_above_hull - decomposition_enthalpy.clip(0) > 0.1"
)
# ax.figure.savefig(f"{module_dir}/mp-e-above-hull-vs-decomp-enth.webp", dpi=300)
