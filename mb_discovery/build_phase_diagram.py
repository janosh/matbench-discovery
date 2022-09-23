# %%
import gzip
import json
import os
import pickle
from datetime import datetime

import pandas as pd
import pymatviz
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.ext.matproj import MPRester

from mb_discovery import ROOT
from mb_discovery.compute_formation_energy import (
    get_elemental_ref_entries,
    get_form_energy_per_atom,
)

today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)


# %%
all_mp_computed_structure_entries = MPRester().get_entries("")  # run on 2022-09-16

# save all ComputedStructureEntries to disk
pd.Series(
    {e.entry_id: e for e in all_mp_computed_structure_entries}
).drop_duplicates().to_json(  # mp-15590 appears twice so we drop_duplicates()
    f"{ROOT}/data/{today}-all-mp-entries.json.gz", default_handler=lambda x: x.as_dict()
)


# %%
all_mp_computed_entries = (
    pd.read_json(f"{ROOT}/data/2022-09-16-all-mp-entries.json.gz")
    .set_index("material_id")
    .entry.map(ComputedEntry.from_dict)  # drop the structure, just load ComputedEntry
    .to_dict()
)


print(f"{len(all_mp_computed_entries) = :,}")
# len(all_mp_computed_entries) = 146,323


# %% build phase diagram with MP entries only
ppd_mp = PatchedPhaseDiagram(all_mp_computed_entries)
# prints:
# PatchedPhaseDiagram
#   Covering 44805 Sub-Spaces

# save MP PPD to disk
with gzip.open(f"{module_dir}/{today}-ppd-mp.pkl.gz", "wb") as zip_file:
    pickle.dump(ppd_mp, zip_file)


# %% build phase diagram with both MP entries + WBM entries
df_wbm = pd.read_json(
    f"{ROOT}/data/2022-06-26-wbm-cses-and-initial-structures.json.gz"
).set_index("material_id")

wbm_computed_entries: list[ComputedEntry] = df_wbm.query("n_elements > 1").cse.map(
    ComputedEntry.from_dict
)

wbm_computed_entries = MaterialsProject2020Compatibility().process_entries(
    wbm_computed_entries, verbose=True, clean=True
)

n_skipped = len(df_wbm) - len(wbm_computed_entries)
print(f"{n_skipped:,} ({n_skipped / len(df_wbm):.1%}) entries not processed")


# %% merge MP and WBM entries into a single PatchedPhaseDiagram
mp_wbm_ppd = PatchedPhaseDiagram(
    wbm_computed_entries + all_mp_computed_entries, verbose=True
)


# %% compute terminal reference entries across all MP (can be used to compute MP
# compatible formation energies quickly)
elemental_ref_entries = get_elemental_ref_entries(all_mp_computed_entries)

# save elemental_ref_entries to disk as json
with open(f"{module_dir}/{today}-elemental-ref-entries.json", "w") as file:
    json.dump(elemental_ref_entries, file, default=lambda x: x.as_dict())


# %% load MP elemental reference entries to compute formation energies
mp_elem_refs_path = f"{ROOT}/data/2022-09-19-mp-elemental-reference-entries.json"
mp_reference_entries = (
    pd.read_json(mp_elem_refs_path, typ="series").map(ComputedEntry.from_dict).to_dict()
)


df_mp = pd.read_json(f"{ROOT}/data/2022-08-13-mp-all-energies.json.gz").set_index(
    "material_id"
)


# %%
df_mp["our_mp_e_form"] = [
    get_form_energy_per_atom(all_mp_computed_entries[mp_id], mp_reference_entries)
    for mp_id in df_mp.index
]


# make sure get_form_energy_per_atom() reproduces MP formation energies
ax = pymatviz.density_scatter(
    df_mp["formation_energy_per_atom"], df_mp["our_mp_e_form"]
)
ax.set(
    title="MP Formation Energy Comparison",
    xlabel="MP Formation Energy (eV/atom)",
    ylabel="Our Formation Energy (eV/atom)",
)
ax.figure.savefig(f"{ROOT}/tmp/{today}-mp-formation-energy-comparison.png", dpi=300)
