# %%
import gzip
import json
import warnings

import pandas as pd
from pymatgen.entries.compatibility import (
    MaterialsProject2020Compatibility,
    MaterialsProjectCompatibility,
)
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from tqdm import tqdm

from matbench_discovery import ROOT, today
from matbench_discovery.data import df_wbm
from matbench_discovery.energy import get_e_form_per_atom
from matbench_discovery.plots import plt

"""
NOTE MaterialsProject2020Compatibility takes structural information into account when
correcting energies (only applies to certain oxides and sulfides). Always use
ComputedStructureEntry, not ComputedEntry when applying corrections.
"""


cse_path = f"{ROOT}/data/wbm/2022-10-19-wbm-computed-structure-entries.json.bz2"
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

df_wbm["mp2020_cse_correction_per_atom"] = [
    cse.correction_per_atom for cse in tqdm(cses)
]
df_wbm["mp2020_ce_correction_per_atom"] = [ce.correction_per_atom for ce in tqdm(ces)]


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
ax = plt.gca()
for key, df_anion in df_ce_ne_cse.groupby("anion"):
    ax = df_anion.plot.scatter(
        ax=ax,
        x="mp2020_cse_correction_per_atom",
        y="mp2020_ce_correction_per_atom",
        label=f"{key} ({len(df_anion):,})",
        color=dict(oxide="orange", sulfide="teal").get(key, "blue"),
        title=f"CSE vs CE corrections for ({len(df_ce_ne_cse):,} / {len(df_wbm):,} = "
        f"{len(df_ce_ne_cse) / len(df_wbm):.1%})\n outliers of largest difference",
    )

ax.axline((0, 0), slope=1, color="gray", linestyle="dashed", zorder=-1)

ax.figure.savefig(f"{ROOT}/tmp/{today}-ce-vs-cse-corrections-outliers.pdf")


# %%
ax = plt.gca()
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
ax.figure.savefig(f"{ROOT}/tmp/{today}-ce-vs-cse-e-form-outliers.pdf")


# %% below code resulted in
# https://github.com/materialsproject/pymatgen/issues/2730
wbm_step_2_34803 = (
    df_ce_ne_cse.e_form_per_atom_mp2020_from_cse
    - df_ce_ne_cse.e_form_per_atom_mp2020_from_ce
).idxmax()
idx = df_wbm.index.get_loc(wbm_step_2_34803)
cse_mp2020, cse_legacy = cses[idx].copy(), cses[idx].copy()
ce_mp2020, ce_legacy = ces[idx].copy(), ces[idx].copy()


with gzip.open(f"{ROOT}/tmp/cse-wbm-2-34803.json.zip", "w") as f:
    f.write(cse_mp2020.to_json().encode("utf-8"))

with gzip.open(f"{ROOT}/tmp/cse-wbm-2-34803.json.zip") as f:
    cse = ComputedStructureEntry.from_dict(json.load(f))

cse_mp2020 = cse.copy()
cse_legacy = cse.copy()
ce_mp2020 = ComputedEntry.from_dict(cse.to_dict())
ce_legacy = ce_mp2020.copy()


MaterialsProject2020Compatibility().process_entry(cse_mp2020)
MaterialsProject2020Compatibility().process_entry(ce_mp2020)
MaterialsProjectCompatibility().process_entry(cse_legacy)
MaterialsProjectCompatibility().process_entry(ce_legacy)

print(f"{cse_mp2020.correction=:.4}")
print(f"{ce_mp2020.correction=:.4}")
print(f"{cse_legacy.correction=:.4}")
print(f"{ce_legacy.correction=:.4}")

print(f"{cse_mp2020.energy_adjustments=}\n")
print(f"{ce_mp2020.energy_adjustments=}\n")
print(f"{cse_legacy.energy_adjustments=}\n")
print(f"{ce_legacy.energy_adjustments=}\n")
