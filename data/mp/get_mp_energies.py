"""Download all MP formation and above hull energies on 2023-01-10.

The main purpose of this script is produce the file at DataFiles.mp_energies.path.

Related EDA of MP formation energies:
https://github.com/janosh/pymatviz/blob/-/examples/mp_bimodal_e_form.ipynb
"""

# %%
import os

import pandas as pd
import plotly.express as px
import pymatviz as pmv
from mp_api.client import MPRester
from pymatgen.core import Element, Structure
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import STABILITY_THRESHOLD, today
from matbench_discovery.data import DataFiles
from matbench_discovery.structure import prototype

__author__ = "Janosh Riebesell"
__date__ = "2023-01-10"

module_dir = os.path.dirname(__file__)
potcar_spec_key = "potcar_spec"
fields = {
    Key.mat_id,
    "formula_pretty",
    e_form_key := Key.formation_energy_per_atom,
    e_per_atom_key := "energy_per_atom",
    "symmetry",
    "energy_above_hull",
    e_decomp_enth_key := "decomposition_enthalpy",
    "energy_type",
    "nsites",
    thermo_type_key := "thermo_type",
}


# %%
with MPRester(use_document_model=False) as mpr:
    docs = mpr.thermo.search(fields=fields, thermo_types=["GGA_GGA+U"])

assert fields == set(docs[0]), f"missing fields: {fields - set(docs[0])}"
print(f"{today}: {len(docs)=:,}")
# 2022-08-13: len(docs) = 146,323
# 2023-01-10: len(docs) = 154,718


# %%
df_mp = pd.DataFrame(docs).set_index(Key.mat_id)
df_mp = df_mp.rename(columns={"formula_pretty": Key.formula, "nsites": Key.n_sites})

df_spg = pd.json_normalize(df_mp.pop("symmetry"))[["number", "symbol"]]
df_mp["spacegroup_symbol"] = df_spg.symbol.to_numpy()

df_mp.energy_type.value_counts().plot.pie(autopct="%1.1f%%")
# GGA: 72.2%, GGA+U: 27.8%


# %%
df_mp_cse = pd.read_json(DataFiles.mp_computed_structure_entries.path).set_index(
    Key.mat_id
)

df_mp_cse[Key.structure] = [
    Structure.from_dict(cse[Key.structure])
    for cse in tqdm(df_mp_cse.entry, desc="Hydrating structures")
]
df_mp_cse[f"{Key.protostructure}_moyo"] = [
    prototype.get_protostructure_label(struct)
    for struct in tqdm(df_mp_cse.structure, desc="Calculating proto-structure labels")
]
# make sure symmetry detection succeeded for all structures
assert df_mp_cse[f"{Key.protostructure}_moyo"].str.startswith("invalid").sum() == 0
df_mp[f"{Key.protostructure}_moyo"] = df_mp_cse[f"{Key.protostructure}_moyo"]

spg_nums = df_mp[f"{Key.protostructure}_moyo"].str.split("_").str[2].astype(int)
# make sure all our spacegroup numbers match MP's
assert (spg_nums.sort_index() == df_spg["number"].sort_index()).all()

df_mp.to_csv(DataFiles.mp_energies.path)
# df_mp = pd.read_csv(DataFiles.mp_energies.path, na_filter=False).set_index(Key.mat_id)


# %% reproduce fig. 1b from https://arxiv.org/abs/2001.10591 (as data consistency check)
ax = df_mp.plot.scatter(
    x=e_form_key,
    y=e_decomp_enth_key,
    alpha=0.1,
    xlim=[-5, 1],
    ylim=[-1, 1],
    color=(df_mp[e_decomp_enth_key] > STABILITY_THRESHOLD).map(
        {True: "red", False: "blue"}
    ),
    title=f"{today} - {len(df_mp):,} MP entries",
)

pmv.powerups.annotate_metrics(df_mp[e_form_key], df_mp[e_decomp_enth_key])
# result on 2023-01-10: plots match. no correlation between formation energy and
# decomposition enthalpy. R^2 = -1.571, MAE = 1.604
# pmv.save_fig(ax, f"{module_dir}/mp-decomp-enth-vs-e-form.webp", dpi=300)


# %% scatter plot energy above convex hull vs decomposition enthalpy
# https://berkeleytheory.slack.com/archives/C16RE1TUN/p1673887564955539
mask_above_line = df_mp.energy_above_hull - df_mp[e_decomp_enth_key].clip(0) > 0.1
ax = df_mp.plot.scatter(
    x=e_decomp_enth_key,
    y="energy_above_hull",
    color=mask_above_line.map({True: "red", False: "blue"}),
    hover_data=["index", Key.formula, e_form_key],
)
# most points lie on line y=x for x > 0 and y = 0 for x < 0.
n_above_line = sum(mask_above_line)
ax.set(
    title=f"{n_above_line:,} / {len(df_mp):,} = {n_above_line / len(df_mp):.1%} "
    f"MP materials with\nenergy_above_hull - {e_decomp_enth_key}.clip(0) > 0.1"
)
# pmv.save_fig(ax, f"{module_dir}/mp-e-above-hull-vs-decomp-enth.webp", dpi=300)


# %% pull all single-element entries and filter for lowest energy per atom system for
# each thermo type to use as elemental reference energies
with MPRester(use_document_model=False) as mpr:
    docs = mpr.thermo.search(num_elements=1)


# %%
df_elem_refs = pd.DataFrame(docs).set_index(Key.mat_id, drop=False)
df_elem_refs[Key.atomic_number] = [Element(el).Z for el in df_elem_refs["chemsys"]]
df_elem_refs[potcar_spec_key] = [
    ", ".join(
        potcar["titel"]  # codespell:ignore
        for potcar in next(iter(entry.values())).parameters["potcar_spec"]
    )
    for entry in df_elem_refs["entries"]
]


# groupby element and take min energy per atom
df_elem_refs = df_elem_refs.sort_values(by=e_per_atom_key)
df_elem_refs = (
    df_elem_refs.groupby(["formula_pretty", thermo_type_key]).first().reset_index()
)


# %%
fig = px.scatter(
    df_elem_refs,
    x=Key.atomic_number,
    y=e_per_atom_key,
    color=thermo_type_key,
    hover_name=Key.mat_id,
    hover_data=["formula_pretty", "energy_above_hull", potcar_spec_key],
)
fig.show()
