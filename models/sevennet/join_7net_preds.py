from glob import glob

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import DataFiles, as_dict_handler, df_wbm
from matbench_discovery.energy import get_e_form_per_atom, mp_elemental_ref_energies

e_form_7net_col = "e_form_per_atom_sevennet"
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

df_7net = pd.concat(dfs.values()).round(4)

if len(df_7net) != len(df_wbm):  # make sure there is no missing structure
    raise ValueError("Missing structures in SevenNet results")

df_cse = pd.read_json(DataFiles.wbm_computed_structure_entries.path).set_index(
    Key.mat_id
)
df_cse[Key.cse] = [
    ComputedStructureEntry.from_dict(dct) for dct in tqdm(df_cse[Key.cse])
]

# As SevenNet-0 (11July2024) is trained on 'uncorrected energy' of MPtrj,
# MP formation energy corrections need to be applied


# %% transfer energies and relaxed structures WBM CSEs since MP2020 energy
# corrections applied below are structure-dependent (for oxides and sulfides)
cse: ComputedStructureEntry
for row in tqdm(df_7net.itertuples(), total=len(df_7net), desc="ML energies to CSEs"):
    mat_id, struct_dict, energy, *_ = row
    mlip_struct = Structure.from_dict(struct_dict)
    cse = df_cse.loc[mat_id, Key.cse]
    cse._energy = energy  # cse._energy is the uncorrected energy from MPtrj dataset (or vasp raw)  # noqa: E501, SLF001
    cse._structure = mlip_struct  # noqa: SLF001
    df_7net.loc[mat_id, Key.cse] = cse


orig_len = len(df_7net)
MaterialsProject2020Compatibility().process_entries(  # change in-place
    df_7net[Key.cse], verbose=True, clean=True
)
if len(df_7net) != orig_len:
    raise ValueError("Some structures were removed during energy correction")

df_7net[Key.formula] = df_wbm[Key.formula]

df_7net[e_form_7net_col] = [
    # see https://matbench-discovery.materialsproject.org/data#mp-elemental-reference-energies
    # MP ref energies are the lowest energies found for unary structures of each element
    get_e_form_per_atom(
        dict(energy=cse.energy, composition=formula), mp_elemental_ref_energies
    )
    for formula, cse in tqdm(
        df_7net.set_index(Key.formula)[Key.cse].items(), total=len(df_7net)
    )
]
df_7net = df_7net.round(4)

df_7net.select_dtypes("number").to_csv(f"{out_path}.csv.gz")  # save csv storable
df_7net.reset_index().to_json(f"{out_path}.json.gz", default_handler=as_dict_handler)
