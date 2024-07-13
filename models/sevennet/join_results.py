import sys
from glob import glob

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from tqdm import tqdm

from matbench_discovery.data import DATA_FILES, as_dict_handler, df_wbm
from matbench_discovery.energy import get_e_form_per_atom, mp_elemental_ref_energies
from matbench_discovery.enums import Key

STRUCT_COL = "sevennet_structure"
E_FORM_COL = "e_form_per_atom_sevennet"

results = "./results"
pot_name = "sevennet"
out_path = f"{results}/{pot_name}"
files = sorted(glob(f"{results}/{pot_name}-*.json.gz"))

dfs = {}
for file_path in tqdm(files):
    if file_path in dfs:
        continue
    df_i = pd.read_json(file_path).set_index(Key.mat_id)
    dfs[file_path] = df_i

df_sevenn = pd.concat(dfs.values()).round(4)

assert len(df_sevenn) == len(df_wbm)  # make sure there is no missing structure

df_cse = pd.read_json(DATA_FILES.wbm_computed_structure_entries).set_index(Key.mat_id)
df_cse[Key.cse] = [
    ComputedStructureEntry.from_dict(dct) for dct in tqdm(df_cse[Key.cse])
]

# As SevenNet-0 (11July2024) is trained on 'uncorrected energy' of MPTrj,
# energies should be corrected basd on MP

# %% transfer energies and relaxed structures WBM CSEs since MP2020 energy
# corrections applied below are structure-dependent (for oxides and sulfides)
cse: ComputedStructureEntry  # pymatgen class that has both structure and computed energy
for row in tqdm(df_sevenn.itertuples(), total=len(df_sevenn), desc="ML energies to CSEs"):
    mat_id, struct_dict, energy, *_ = row
    mlip_struct = Structure.from_dict(struct_dict)
    df_sevenn.at[mat_id, STRUCT_COL] = mlip_struct  # noqa: PD008
    cse = df_cse.loc[mat_id, Key.cse]
    cse._energy = energy  # cse._energy is the uncorrected energy from MPtraj dataset (or vasp raw)
    cse._structure = mlip_struct  # noqa: SLF001
    df_sevenn.loc[mat_id, Key.cse] = cse


orig_len = len(df_sevenn)
MaterialsProject2020Compatibility().process_entries(  # change in-place
    df_sevenn[Key.cse], verbose=True, clean=True
)
assert len(df_sevenn) == orig_len

############# df_sevenn is now have corrected energy #############

df_sevenn[Key.formula] = df_wbm[Key.formula]  # pd automatically match rows by material_id

df_sevenn[E_FORM_COL] = [
    # see https://matbench-discovery.materialsproject.org/data#mp-elemental-reference-energies
    # MP ref energies are 'the lowest energies found for unary structures of each element'
    get_e_form_per_atom(dict(energy=cse.energy, composition=formula), mp_elemental_ref_energies)
    for formula, cse in tqdm(
        df_sevenn.set_index(Key.formula)[Key.cse].items(), total=len(df_sevenn)
    )
]
df_wbm[E_FORM_COL] = df_sevenn[E_FORM_COL]

df_sevenn = df_sevenn.round(4)

df_sevenn.select_dtypes("number").to_csv(f"{out_path}.csv.gz")  # save csv storable
df_sevenn.reset_index().to_json(f"{out_path}.json.gz", default_handler=as_dict_handler)
