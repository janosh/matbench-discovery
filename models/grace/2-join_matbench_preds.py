"""Join and process GRACE model predictions for matbench-discovery.

This processing script has been copied from the nequip script here:
https://github.com/janosh/matbench-discovery/blob/main/models/nequip/join_matbench_preds.py
And then slightly refactored for GRACE.

Uses matbench-discovery commit ID 012ccfe, k_srme commit ID 0269a946,
pymatviz v0.15.1.
"""

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

model_name = "grace_2l_oam_l"
e_form_potential_col = f"e_form_per_atom_{model_name}"
results_dir = "./results"

out_path = f"{results_dir}/{model_name}"
files = sorted(glob(f"{results_dir}/{model_name}/*.json.gz"))

dfs_loaded = {}
for file_path in tqdm(files, desc="Loading results"):
    if file_path in dfs_loaded:
        continue
    df_result = pd.read_json(file_path, orient="records", lines=True).set_index(
        Key.mat_id
    )
    dfs_loaded[file_path] = df_result

df_model = pd.concat(dfs_loaded.values()).round(4)

if len(df_model) != len(df_wbm):  # make sure there is no missing structure
    warnings.warn(
        f"Some missing structures in results, {len(df_model)} in df_model, "
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
computed_struct_entry: ComputedStructureEntry
for row in tqdm(
    df_model.itertuples(),
    total=len(df_model),
    desc="Generating ML-predicted ComputedStructureEntrys",
):
    mat_id, struct_dict, energy, *_rest = row
    mlip_struct = Structure.from_dict(struct_dict)
    computed_struct_entry = df_cse.loc[mat_id, Key.computed_structure_entry]
    # computed_struct_entry._energy is the uncorrected energy from MPtrj dataset
    computed_struct_entry._energy = energy  # noqa: SLF001
    computed_struct_entry._structure = mlip_struct  # noqa: SLF001
    df_model.loc[mat_id, Key.computed_structure_entry] = computed_struct_entry


orig_len = len(df_model)
print("Applying MP 2020 Compatibility corrections")
df_model[Key.computed_structure_entry] = (
    MaterialsProject2020Compatibility().process_entries(  # change in-place
        df_model[Key.computed_structure_entry],
        verbose=True,
        clean=True,
        inplace=False,
        n_workers=8,  # faster processing
    )
)
if len(df_model) != orig_len:
    raise ValueError("Some structures were removed during energy correction")

df_model[Key.formula] = df_wbm[Key.formula]

df_model[e_form_potential_col] = [
    # see https://matbench-discovery.materialsproject.org/data#mp-elemental-reference-energies
    # MP ref energies are the lowest energies found for unary structures of each element
    calc_energy_from_e_refs(
        dict(energy=entry.energy, composition=formula), mp_elemental_ref_energies
    )
    for formula, entry in tqdm(
        df_model.set_index(Key.formula)[Key.computed_structure_entry].items(),
        total=len(df_model),
        desc="Computing MP-corrected formation energies",
    )
]
df_model = df_model.round(4)

csv_out_path = f"{out_path}.csv.gz"
print(f"Saving to file {csv_out_path}")
df_model.select_dtypes("number").to_csv(csv_out_path)

json_out_path = f"{out_path}.json.gz"
print(f"Saving to file {json_out_path}")
df_model.reset_index().to_json(json_out_path, default_handler=as_dict_handler)
