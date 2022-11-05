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
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.ext.matproj import MPRester

from matbench_discovery import ROOT
from matbench_discovery.energy import get_e_form_per_atom, get_elemental_ref_entries

today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)


# %%
all_mp_computed_structure_entries = MPRester().get_entries("")  # run on 2022-09-16

# save all ComputedStructureEntries to disk
pd.Series(
    {e.entry_id: e for e in all_mp_computed_structure_entries}
).drop_duplicates().to_json(  # mp-15590 appears twice so we drop_duplicates()
    f"{module_dir}/{today}-mp-computed-structure-entries.json.gz",
    default_handler=lambda x: x.as_dict(),
)


# %%
data_path = f"{module_dir}/2022-09-16-mp-computed-structure-entries.json.gz"
df = pd.read_json(data_path).set_index("material_id")
# drop the structure, just load ComputedEntry
mp_computed_entries = df.entry.map(ComputedEntry.from_dict).to_dict()

print(f"{len(mp_computed_entries) = :,}")
# len(mp_computed_entries) = 146,323


# %% build phase diagram with MP entries only
ppd_mp = PatchedPhaseDiagram(mp_computed_entries)
# prints:
# PatchedPhaseDiagram covering 44805 sub-spaces

# save MP PPD to disk
with gzip.open(f"{module_dir}/{today}-ppd-mp.pkl.gz", "wb") as zip_file:
    pickle.dump(ppd_mp, zip_file)


# %% build phase diagram with both MP entries + WBM entries
df_wbm = pd.read_json(
    f"{ROOT}/data/wbm/2022-10-19-wbm-cses+init-structs.json.bz2"
).set_index("material_id")

# using ComputedStructureEntry vs ComputedEntry here is important as CSEs receive
# more accurate energy corrections that take into account peroxide/superoxide nature
# of materials (and same for sulfides) based on atomic distances in the structure
wbm_computed_entries: list[ComputedStructureEntry] = df_wbm.query(
    "n_elements > 1"
).cse.map(ComputedStructureEntry.from_dict)

wbm_computed_entries = MaterialsProject2020Compatibility().process_entries(
    wbm_computed_entries, verbose=True, clean=True
)

n_skipped = len(df_wbm) - len(wbm_computed_entries)
assert n_skipped == 0
print(f"{n_skipped:,} ({n_skipped / len(df_wbm):.1%}) entries not processed")


# %% merge MP and WBM entries into a single PatchedPhaseDiagram
mp_wbm_ppd = PatchedPhaseDiagram(
    wbm_computed_entries + mp_computed_entries, verbose=True
)

# save MP+WBM PPD to disk (not run)
with gzip.open(f"{module_dir}/{today}-ppd-mp.pkl.gz", "wb") as zip_file:
    pickle.dump(mp_wbm_ppd, zip_file)


# %% compute terminal reference entries across all MP (can be used to compute MP
# compatible formation energies quickly)
elemental_ref_entries = get_elemental_ref_entries(mp_computed_entries)

# save elemental_ref_entries to disk as json
with open(f"{ROOT}/data/mp/{today}-mp-elemental-reference-entries.json", "w") as file:
    json.dump(elemental_ref_entries, file, default=lambda x: x.as_dict())


df_mp = pd.read_json(f"{ROOT}/data/mp/2022-08-13-mp-energies.json.gz").set_index(
    "material_id"
)


# %%
df_mp["our_mp_e_form"] = [
    get_e_form_per_atom(mp_computed_entries[mp_id]) for mp_id in df_mp.index
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
ax.figure.savefig(f"{ROOT}/tmp/{today}-our-vs-mp-formation-energies.png", dpi=300)
