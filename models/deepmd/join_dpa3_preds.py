"""Join DPA3 predictions from multiple files and compute formation energies.

This script combines DPA3 model predictions from multiple JSON files, applies MP2020
energy corrections, and computes formation energies using MP elemental references.
"""

import os
from glob import glob

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import df_wbm
from matbench_discovery.energy import get_e_form_per_atom, mp_elemental_ref_energies
from matbench_discovery.enums import DataFiles

e_form_dp_col = "e_form_per_atom_dp"
results = "./results"
model_name = "dpa3"
module_dir = os.path.dirname(__file__)
out_path = f"{module_dir}/{model_name}"
files = sorted(glob(f"{results}/{model_name}-*.json.gz"))

dfs = {}
for file_path in tqdm(files):
    if file_path in dfs:
        continue
    df_i = pd.read_json(file_path).set_index(Key.mat_id)
    dfs[file_path] = df_i

df_dpa3 = pd.concat(dfs.values()).round(4)

if len(df_dpa3) != len(df_wbm):  # make sure there is no missing structure
    raise ValueError("Missing structures in DPA3 results")

df_cse = pd.read_json(DataFiles.wbm_computed_structure_entries.path).set_index(
    Key.mat_id
)
df_cse[Key.computed_structure_entry] = [
    ComputedStructureEntry.from_dict(dct)
    for dct in tqdm(df_cse[Key.computed_structure_entry], desc="Hydrate CSEs")
]

# As DPA3 is trained on 'uncorrected energy' of MPtrj,
# MP formation energy corrections need to be applied


# transfer energies and relaxed structures WBM CSEs since MP2020 energy
# corrections applied below are structure-dependent (for oxides and sulfides)
cse: ComputedStructureEntry
for row in tqdm(df_dpa3.itertuples(), total=len(df_dpa3), desc="ML energies to CSEs"):
    mat_id, struct_dict, energy, *_ = row
    mlip_struct = Structure.from_dict(struct_dict)
    cse = df_cse.loc[mat_id, Key.computed_structure_entry]
    cse._energy = energy  # cse._energy is the uncorrected energy from MPtrj dataset (or vasp raw)  # noqa: E501, SLF001
    cse._structure = mlip_struct  # noqa: SLF001
    df_dpa3.loc[mat_id, Key.computed_structure_entry] = cse


orig_len = len(df_dpa3)
MaterialsProject2020Compatibility().process_entries(  # change in-place
    df_dpa3[Key.computed_structure_entry], verbose=True, clean=True
)
if len(df_dpa3) != orig_len:
    raise ValueError("Some structures were removed during energy correction")

df_dpa3[Key.formula] = df_wbm[Key.formula]

df_dpa3[e_form_dp_col] = [
    # see https://matbench-discovery.materialsproject.org/data#mp-elemental-reference-energies
    # MP ref energies are the lowest energies found for unary structures of each element
    get_e_form_per_atom(
        dict(energy=cse.energy, composition=formula), mp_elemental_ref_energies
    )
    for formula, cse in tqdm(
        df_dpa3.set_index(Key.formula)[Key.computed_structure_entry].items(),
        total=len(df_dpa3),
    )
]
df_dpa3 = df_dpa3.round(4)

df_dpa3.select_dtypes("number").to_csv(f"{out_path}.csv.gz")  # save CSV storable
# df_dpa3.reset_index().to_json(f"{out_path}.json.gz", default_handler=as_dict_handler)
