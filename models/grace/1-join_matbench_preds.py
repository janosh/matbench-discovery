"""
This processing script has been copied from the nequip script here:
https://github.com/janosh/matbench-discovery/blob/main/models/nequip/join_matbench_preds.py
And then slightly refactored for GRACE
"""

# uses matbench-discovery matbench-discovery commit ID 012ccfe,
# k_srme commit ID 0269a946, pymatviz v0.15.1
import warnings
from glob import glob

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import DataFiles, as_dict_handler, df_wbm
from matbench_discovery.energy import calc_energy_from_e_refs, mp_elemental_ref_energies

pot_name = "grace_2l_oam_l"
pot_name = pot_name.lower()
e_form_potential_col = f"e_form_per_atom_{pot_name}"
results = "./results"

out_path = f"{results}/{pot_name}"
files = sorted(glob(f"{results}/{pot_name}/*.json.gz"))

dfs = {}
for file_path in tqdm(files, desc="Loading results"):
    if file_path in dfs:
        continue
    df_i = pd.read_json(file_path, orient="records", lines=True).set_index(Key.mat_id)
    dfs[file_path] = df_i

df_potential = pd.concat(dfs.values()).round(4)

if len(df_potential) != len(df_wbm):  # make sure there is no missing structure
    warnings.warn(
        f"Some missing structures in results, {len(df_potential)} in df_potential, "
        f"{len(df_wbm)} in df_wbm, likely due to some crashed relaxations",
        stacklevel=2,
    )
    # raise ValueError("Missing structures in results")

print("Loading reference WBM dataset")
df_cse = pd.read_json(
    DataFiles.wbm_computed_structure_entries.path, orient="records", lines=True
).set_index(Key.mat_id)
df_cse[Key.computed_structure_entry] = [
    ComputedStructureEntry.from_dict(dct)
    for dct in tqdm(
        df_cse[Key.computed_structure_entry],
        desc="Generating WBM reference ComputedStructureEntrys",
    )
]

# trained on 'uncorrected energy' of MPtrj,
# MP formation energy corrections need to be applied


# %% transfer energies and relaxed structures WBM CSEs since MP2020 energy
# corrections applied below are structure-dependent (for oxides and sulfides)
cse: ComputedStructureEntry
for row in tqdm(
    df_potential.itertuples(),
    total=len(df_potential),
    desc="Generating ML-predicted ComputedStructureEntrys",
):
    mat_id, struct_dict, energy, *_ = row
    mlip_struct = Structure.from_dict(struct_dict)
    cse = df_cse.loc[mat_id, Key.computed_structure_entry]
    cse._energy = energy  # cse._energy is the uncorrected energy from MPtrj dataset (or vasp raw)  # noqa: E501, SLF001
    cse._structure = mlip_struct  # noqa: SLF001
    df_potential.loc[mat_id, Key.computed_structure_entry] = cse


orig_len = len(df_potential)
print("Applying MP 2020 Compatibility corrections")
df_potential[Key.computed_structure_entry] = (
    MaterialsProject2020Compatibility().process_entries(  # change in-place
        df_potential[Key.computed_structure_entry],
        verbose=True,
        clean=True,
        inplace=False,
        n_workers=8,  # faster processing
    )
)
if len(df_potential) != orig_len:
    raise ValueError("Some structures were removed during energy correction")

df_potential[Key.formula] = df_wbm[Key.formula]

df_potential[e_form_potential_col] = [
    # see https://matbench-discovery.materialsproject.org/data#mp-elemental-reference-energies
    # MP ref energies are the lowest energies found for unary structures of each element
    calc_energy_from_e_refs(
        dict(energy=cse.energy, composition=formula), mp_elemental_ref_energies
    )
    for formula, cse in tqdm(
        df_potential.set_index(Key.formula)[Key.computed_structure_entry].items(),
        total=len(df_potential),
        desc="Getting formation energies",
    )
]
df_potential = df_potential.round(4)

print(f"Saving to file {out_path}.csv.gz")
df_potential.select_dtypes("number").to_csv(f"{out_path}.csv.gz")  # save csv storable
# we can also save the full computed structure entries to file if desired:
print(f"Saving to file {out_path}.json.gz")
df_potential.reset_index().to_json(
    f"{out_path}.json.gz", default_handler=as_dict_handler
)
