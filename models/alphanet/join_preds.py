from glob import glob

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import as_dict_handler, df_wbm
from matbench_discovery.energy import get_e_form_per_atom, mp_elemental_ref_energies
from matbench_discovery.enums import DataFiles

e_form_Anet_col = "e_form_per_atom_alphanet"  # noqa: N816
results = "./res_relax"
pot_name = "alphanet"
out_path = f"{results}/{pot_name}"
files = sorted(glob(f"{results}/{pot_name}-*.json.gz"))

dfs = {}
for file_path in tqdm(files):
    if file_path in dfs:
        continue
    df_i = pd.read_json(file_path).set_index(Key.mat_id)
    dfs[file_path] = df_i

df_Anet = pd.concat(dfs.values()).round(4)  # noqa: N816

if len(df_Anet) != len(df_wbm):  # make sure there is no missing structure
    raise ValueError("Missing structures in SevenNet results")

wbm_cse_path = DataFiles.wbm_computed_structure_entries.path
df_cse = pd.read_json(wbm_cse_path).set_index(Key.mat_id)

df_cse[Key.computed_structure_entry] = [
    ComputedStructureEntry.from_dict(dct)
    for dct in tqdm(df_cse[Key.computed_structure_entry], desc="Hydrate CSEs")
]


# %% transfer energies and relaxed structures WBM CSEs since MP2020 energy
# corrections applied below are structure-dependent (for oxides and sulfides)
cse: ComputedStructureEntry
for row in tqdm(df_Anet.itertuples(), total=len(df_Anet), desc="ML energies to CSEs"):
    mat_id, struct_dict, energy, *_ = row
    mlip_struct = Structure.from_dict(struct_dict)
    cse = df_cse.loc[mat_id, Key.computed_structure_entry]
    cse._energy = energy  # noqa: SLF001
    cse._structure = mlip_struct  # noqa: SLF001
    df_Anet.loc[mat_id, Key.computed_structure_entry] = cse


orig_len = len(df_Anet)
MaterialsProject2020Compatibility().process_entries(
    df_Anet[Key.computed_structure_entry], verbose=True, clean=True
)
if len(df_Anet) != orig_len:
    raise ValueError("Some structures were removed during energy correction")

df_Anet[Key.formula] = df_wbm[Key.formula]

df_Anet[e_form_Anet_col] = [
    # see https://matbench-discovery.materialsproject.org/data#mp-elemental-reference-energies
    # MP ref energies are the lowest energies found for unary structures of each element
    get_e_form_per_atom(
        dict(energy=cse.energy, composition=formula), mp_elemental_ref_energies
    )
    for formula, cse in tqdm(
        df_Anet.set_index(Key.formula)[Key.computed_structure_entry].items(),
        total=len(df_Anet),
    )
]
df_Anet = df_Anet.round(4)  # noqa: N816

df_Anet.select_dtypes("number").to_csv(f"{out_path}.csv.gz")  # save csv storable
df_Anet.reset_index().to_json(f"{out_path}.json.gz", default_handler=as_dict_handler)
