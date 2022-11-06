import warnings
from datetime import datetime

import pandas as pd
from pymatgen.entries.compatibility import (
    MaterialsProject2020Compatibility,
    MaterialsProjectCompatibility,
)
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from tqdm import tqdm

from matbench_discovery import ROOT
from matbench_discovery.energy import get_e_form_per_atom
from matbench_discovery.plot_scripts import df_wbm

"""
NOTE MaterialsProject2020Compatibility takes structural information into account when
correcting energies (only applies to certain oxides and sulfides). Always use
ComputedStructureEntry, not ComputedEntry when applying corrections.
"""

today = f"{datetime.now():%Y-%m-%d}"

cse_path = f"{ROOT}/data/wbm/2022-10-19-wbm-cses.json.bz2"
df_cse = pd.read_json(cse_path).set_index("material_id")

cses = [
    ComputedStructureEntry.from_dict(x) for x in tqdm(df_cse.computed_structure_entry)
]

ces = [ComputedEntry.from_dict(x) for x in tqdm(df_cse.computed_structure_entry)]


warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")


# %%
out = MaterialsProject2020Compatibility().process_entries(cses, verbose=True)
assert len(out) == len(df_cse)
out = MaterialsProject2020Compatibility().process_entries(ces, verbose=True)
assert len(out) == len(df_cse)

df_wbm["e_form_per_atom_mp2020_from_ce"] = [
    get_e_form_per_atom(entry) for entry in tqdm(ces)
]
df_wbm["e_form_per_atom_mp2020_from_cse"] = [
    get_e_form_per_atom(entry) for entry in tqdm(cses)
]

df_wbm["mp2020_cse_correction"] = [cse.correction for cse in tqdm(cses)]
df_wbm["mp2020_ce_correction"] = [ce.correction for ce in tqdm(ces)]


# %%
out = MaterialsProjectCompatibility().process_entries(cses, verbose=True)
assert len(out) == len(df_cse)
out = MaterialsProjectCompatibility().process_entries(ces, verbose=True)
assert len(out) == len(df_cse)

df_wbm["e_form_per_atom_legacy_from_ce"] = [
    get_e_form_per_atom(entry) for entry in tqdm(ces)
]
df_wbm["e_form_per_atom_legacy_from_cse"] = [
    get_e_form_per_atom(entry) for entry in tqdm(cses)
]
df_wbm["legacy_cse_correction"] = [cse.correction for cse in tqdm(cses)]
df_wbm["legacy_ce_correction"] = [ce.correction for ce in tqdm(ces)]


# %%
df_wbm["chem_sys"] = df_wbm.formula.str.replace("[0-9]+", "", regex=True).str.split()
df_wbm["anion"] = None
df_wbm["anion"][df_wbm.chem_sys.astype(str).str.contains("'O'")] = "oxide"
df_wbm["anion"][df_wbm.chem_sys.astype(str).str.contains("'S'")] = "sulfide"

assert dict(df_wbm.anion.value_counts()) == {"oxide": 26984, "sulfide": 10606}

df_ce_ne_cse = df_wbm.query(
    "abs(e_form_per_atom_mp2020_from_cse - e_form_per_atom_mp2020_from_ce) > 1e-4"
)


# %%
for key, df_anion in df_ce_ne_cse.groupby("anion"):
    ax = df_anion.plot.scatter(
        ax=locals().get("ax"),
        x="mp2020_cse_correction",
        y="mp2020_ce_correction",
        label=f"{key} ({len(df_anion):,})",
        color=dict(oxide="orange", sulfide="teal").get(key, "blue"),
        title=f"Outliers in formation energy from CSE vs CE ({len(df_ce_ne_cse):,}"
        f" / {len(df_wbm):,} = {len(df_ce_ne_cse) / len(df_wbm):.1%})",
    )

ax.axline((0, 0), slope=1, color="gray", linestyle="dashed", zorder=-1)


# %%
for key, df_anion in df_ce_ne_cse.groupby("anion"):
    ax = df_anion.plot.scatter(
        ax=locals().get("ax"),
        x="e_form_per_atom_mp2020_from_cse",
        y="e_form_per_atom_mp2020_from_ce",
        label=f"{key} ({len(df_anion):,})",
        color=dict(oxide="orange", sulfide="teal").get(key, "blue"),
        title=f"Outliers in formation energy from CSE vs CE ({len(df_ce_ne_cse):,}"
        f" / {len(df_wbm):,} = {len(df_ce_ne_cse) / len(df_wbm):.1%})",
    )

ax.axline((0, 0), slope=1, color="gray", linestyle="dashed", zorder=-1)

# insight: all materials for which ComputedEntry and ComputedStructureEntry give
# different formation energies are oxides or sulfides for which MP 2020 compat takes
# into account structural information to make more accurate corrections.
# ax.figure.savefig(f"{ROOT}/tmp/{today}-ce-vs-cse-outliers.pdf")
