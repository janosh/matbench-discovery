"""Download all MP formation and above hull energies on 2023-01-10.

Related EDA of MP formation energies:
https://github.com/janosh/pymatviz/blob/-/examples/mp_bimodal_e_form.ipynb
"""

# %%
import os

import pandas as pd
import plotly.express as px
import pymatviz as pmv
from aviary.wren.utils import get_protostructure_label_from_spglib
from mp_api.client import MPRester
from pymatgen.core import Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import STABILITY_THRESHOLD, today
from matbench_discovery.data import DataFiles

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
print(f"{today}: {len(docs)=:,}")
# 2022-08-13: len(docs) = 146,323
# 2023-01-10: len(docs) = 154,718
# 2024-08-21: len(docs)=153,243


# %%
df_mp = pd.DataFrame(docs).set_index(Key.mat_id)
df_mp = df_mp.rename(columns={"formula_pretty": Key.formula, "nsites": Key.n_sites})

df_spg = pd.json_normalize(df_mp.pop("symmetry"))[["number", "symbol"]]
df_spg.index = df_mp.index
df_mp["spacegroup_symbol"] = df_spg.symbol.to_numpy()

px.pie(df_mp.energy_type.value_counts().reset_index(), values="count")
# 2023-01-10: GGA: 72.2%, GGA+U: 27.8%
# 2024-08-21: GGA: 71.8%, GGA+U: 28.2%


# %%
df_cse = pd.read_json(DataFiles.mp_computed_structure_entries.path).set_index(
    Key.mat_id
)


# %%
df_cse[Key.structure] = [
    Structure.from_dict(cse[Key.structure]) for cse in tqdm(df_cse.entry)
]


# %%
df_cse[Key.wyckoff] = [
    get_protostructure_label_from_spglib(struct, raise_errors=True)
    for struct in tqdm(df_cse.structure)
]


# %%
# make sure symmetry detection succeeded for all structures
assert df_cse[Key.wyckoff].str.startswith("invalid").sum() == 0
df_mp[Key.wyckoff] = df_cse[Key.wyckoff]

assert df_mp[Key.wyckoff].isna().sum() == 1048
# 2024-08-21: 1,048 structures in the reference MP dataset are now not in MP


# %%
spg_nums = df_mp[Key.wyckoff].str.split("_").str[2].dropna().astype(int)
# make sure all our spacegroup numbers match MP's
assert (
    spg_nums.sort_index() == df_spg.loc[spg_nums.index, "number"].sort_index()
).all()

# df_mp.to_csv(DataFiles.mp_energies.path)
# df = pd.read_csv(DataFiles.mp_energies.path, na_filter=False).set_index(Key.mat_id)


# %% reproduce fig. 1b from https://arxiv.org/abs/2001.10591 (as data consistency check)
ax = df_mp.plot.scatter(
    x=Key.form_energy,
    y="decomposition_enthalpy",
    alpha=0.1,
    xlim=[-5, 1],
    ylim=[-1, 1],
    color=(df_mp.decomposition_enthalpy > STABILITY_THRESHOLD).map(
        {True: "red", False: "blue"}
    ),
    title=f"{today} - {len(df_mp):,} MP entries",
)

pmv.powerups.annotate_metrics(
    df_mp.formation_energy_per_atom,
    df_mp.decomposition_enthalpy,
    fig=ax,
    fmt=".3f",
)
# result on 2023-01-10: plots match. no correlation between formation energy and
# decomposition enthalpy. R^2 = -1.571, MAE = 1.604
# 2024-08-21: R^2 = -1.582, MAE = 1.611
# pmv.save_fig(ax, f"{module_dir}/mp-decomp-enth-vs-e-form.webp", dpi=300)


# %% scatter plot energy above convex hull vs decomposition enthalpy
# https://berkeleytheory.slack.com/archives/C16RE1TUN/p1673887564955539
mask_above_line = df_mp.energy_above_hull - df_mp.decomposition_enthalpy.clip(0) > 0.1
# most points lie on line y=x for x > 0 and y = 0 for x < 0.
n_above_line = sum(mask_above_line)
fig = px.scatter(
    df_mp.reset_index(),
    x="decomposition_enthalpy",
    y="energy_above_hull",
    color=mask_above_line.map({True: "red", False: "blue"}),
    hover_data=[Key.mat_id, Key.formula, Key.form_energy],
    title=(
        f"{n_above_line:,} / {len(df_mp):,} = {n_above_line / len(df_mp):.1%} "
        "MP materials with<br>energy_above_hull - decomposition_enthalpy.clip(0) > 0.1"
    ),
)
# Update layout for title visibility
fig.update_layout(
    margin=dict(t=120),
    autosize=True,
)
fig.show()
# pmv.save_fig(ax, f"{module_dir}/mp-e-above-hull-vs-decomp-enth.webp", dpi=300)
