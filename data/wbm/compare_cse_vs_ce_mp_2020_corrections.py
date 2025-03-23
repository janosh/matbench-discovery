"""MaterialsProject2020Compatibility takes structural information into account when
correcting energies (for oxides and sulfides only). Always pass
ComputedStructureEntry, not ComputedEntry when calling process_entries.
"""

# %%
import gzip
import json

import pandas as pd
import pymatviz as pmv
from pymatgen.entries.compatibility import (
    MaterialsProject2020Compatibility,
    MaterialsProjectCompatibility,
)
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import ROOT, today
from matbench_discovery.data import df_wbm
from matbench_discovery.energy import get_e_form_per_atom
from matbench_discovery.enums import DataFiles

wbm_cse_path = DataFiles.wbm_computed_structure_entries.path
df_wbm_cse = pd.read_json(wbm_cse_path, lines=True).set_index(Key.mat_id)

cses = [
    ComputedStructureEntry.from_dict(dct)
    for dct in tqdm(
        df_wbm_cse[Key.computed_structure_entry],
        desc="Loading ComputedStructureEntries",
    )
]

ces = [
    ComputedEntry.from_dict(dct)
    for dct in tqdm(
        df_wbm_cse[Key.computed_structure_entry], desc="Loading ComputedEntries"
    )
]


# %%
processed = MaterialsProject2020Compatibility().process_entries(cses, verbose=True)
assert len(processed) == len(df_wbm_cse)
processed = MaterialsProject2020Compatibility().process_entries(ces, verbose=True)
assert len(processed) == len(df_wbm_cse)

df_wbm["e_form_per_atom_mp2020_from_ce"] = [
    get_e_form_per_atom(entry)
    for entry in tqdm(ces, desc="Calculating formation energies from ComputedEntries")
]
df_wbm["e_form_per_atom_mp2020_from_cse"] = [
    get_e_form_per_atom(entry)
    for entry in tqdm(
        cses, desc="Calculating formation energies from ComputedStructureEntries"
    )
]

df_wbm["mp2020_cse_correction_per_atom"] = [
    cse.correction_per_atom for cse in tqdm(cses)
]
df_wbm["mp2020_ce_correction_per_atom"] = [ce.correction_per_atom for ce in tqdm(ces)]


# %%
processed = MaterialsProjectCompatibility().process_entries(cses, verbose=True)
assert len(processed) == len(df_wbm_cse)
processed = MaterialsProjectCompatibility().process_entries(ces, verbose=True)
assert len(processed) == len(df_wbm_cse)

df_wbm["e_form_per_atom_legacy_from_ce"] = [
    get_e_form_per_atom(entry) for entry in tqdm(ces)
]
df_wbm["e_form_per_atom_legacy_from_cse"] = [
    get_e_form_per_atom(entry) for entry in tqdm(cses)
]
df_wbm["legacy_cse_correction"] = [cse.correction for cse in tqdm(cses)]
df_wbm["legacy_ce_correction"] = [ce.correction for ce in tqdm(ces)]


# %%
df_wbm[Key.chem_sys] = (
    df_wbm[Key.formula].str.replace("[0-9]+", "", regex=True).str.split()
)
df_wbm["anion"] = None
df_wbm["anion"][df_wbm[Key.chem_sys].astype(str).str.contains("'O'")] = "oxide"
df_wbm["anion"][df_wbm[Key.chem_sys].astype(str).str.contains("'S'")] = "sulfide"

anion_counts = dict(df_wbm.anion.value_counts())
assert anion_counts == {"oxide": 26_984, "sulfide": 10_596}, f"{anion_counts=}"

df_ce_ne_cse = df_wbm.query(
    "abs(e_form_per_atom_mp2020_from_cse - e_form_per_atom_mp2020_from_ce) > 1e-4"
)


# %%
for x_col, y_col, title in (
    ("mp2020_cse_correction_per_atom", "mp2020_ce_correction_per_atom", "correction"),
    ("e_form_per_atom_mp2020_from_cse", "e_form_per_atom_mp2020_from_ce", "e-form"),
):
    fig = df_ce_ne_cse.plot.scatter(
        x=x_col,
        y=y_col,
        color="anion",
        color_discrete_map={"oxide": "orange", "sulfide": "teal"},
        hover_data=[Key.formula],
        backend="plotly",
    )
    title = f"CSE vs CE {title}<br>({len(df_ce_ne_cse):,} / {len(df_wbm):,} "
    title += f"= {len(df_ce_ne_cse) / len(df_wbm):.1%})"
    fig.layout.title.update(text=title, x=0.5, font=dict(size=16))
    fig.layout.margin.t = 40
    fig.layout.legend.update(x=0, title=None)

    pmv.powerups.add_identity_line(fig)

    # Update legend labels with counts
    for trace in fig.data:
        anion = trace.name
        count = len(df_ce_ne_cse[df_ce_ne_cse.anion == anion])
        trace.name = f"{anion} ({count:,})"

    # insight: all materials for which ComputedEntry and ComputedStructureEntry give
    # different formation energies are oxides or sulfides for which MP 2020 compat takes
    # into account structural information to make more accurate corrections.
    pmv.save_fig(fig, f"{ROOT}/tmp/{today}-ce-vs-cse-{title}-outliers.pdf")
    fig.show()


# %% below code resulted in
# https://github.com/materialsproject/pymatgen/issues/2730
wbm_step_2_34803 = (
    df_ce_ne_cse.e_form_per_atom_mp2020_from_cse
    - df_ce_ne_cse.e_form_per_atom_mp2020_from_ce
).idxmax()
idx = df_wbm.index.get_loc(wbm_step_2_34803)
cse_mp2020, cse_legacy = cses[idx].copy(), cses[idx].copy()
ce_mp2020, ce_legacy = ces[idx].copy(), ces[idx].copy()


with gzip.open(f"{ROOT}/tmp/cse-wbm-2-34803.json.zip", mode="w") as file:
    file.write(cse_mp2020.to_json().encode("utf-8"))

with gzip.open(f"{ROOT}/tmp/cse-wbm-2-34803.json.zip") as file:
    cse = ComputedStructureEntry.from_dict(json.load(file))

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
