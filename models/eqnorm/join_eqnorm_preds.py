import warnings
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

e_form_eqnorm_col = "e_form_per_atom_eqnorm"
results = "./results"
pot_name = "eqnorm"
out_path = f"{results}/{pot_name}"
files = sorted(glob(f"{results}/{pot_name}-*.json.gz"))

dfs = {}
for file_path in tqdm(files):
    if file_path in dfs:
        continue
    df_i = pd.read_json(file_path, lines=True).set_index(Key.mat_id)
    dfs[file_path] = df_i

df_eqnorm = pd.concat(dfs.values()).round(4)

if len(df_eqnorm) != len(df_wbm):  # make sure there is no missing structure
    warnings.warn("Missing structures in eqnorm results", stacklevel=2)

wbm_cse_path = DataFiles.wbm_computed_structure_entries.path
df_wbm_cse = pd.read_json(wbm_cse_path, lines=True).set_index(Key.mat_id)

df_wbm_cse[Key.computed_structure_entry] = [
    ComputedStructureEntry.from_dict(dct)
    for dct in tqdm(df_wbm_cse[Key.computed_structure_entry], desc="Hydrate CSEs")
]

# %% transfer energies and relaxed structures WBM CSEs since MP2020 energy
# corrections applied below are structure-dependent (for oxides and sulfides)
cse: ComputedStructureEntry
for row in tqdm(
    df_eqnorm.itertuples(), total=len(df_eqnorm), desc="ML energies to CSEs"
):
    mat_id, struct_dict, energy, *_ = row
    mlip_struct = Structure.from_dict(struct_dict)
    cse = df_wbm_cse.loc[mat_id, Key.computed_structure_entry]
    cse._energy = energy  # cse._energy is the uncorrected energy from MPtrj dataset (or vasp raw)  # noqa: E501, SLF001
    cse._structure = mlip_struct  # noqa: SLF001
    df_eqnorm.loc[mat_id, Key.computed_structure_entry] = cse


orig_len = len(df_eqnorm)
MaterialsProject2020Compatibility().process_entries(  # change in-place
    df_eqnorm[Key.computed_structure_entry], verbose=True, clean=True
)
if len(df_eqnorm) != orig_len:
    warnings.warn("Some structures were removed during energy correction", stacklevel=2)

df_eqnorm[Key.formula] = df_wbm[Key.formula]

df_eqnorm[e_form_eqnorm_col] = [
    # see https://matbench-discovery.materialsproject.org/data#mp-elemental-reference-energies
    # MP ref energies are the lowest energies found for unary structures of each element
    get_e_form_per_atom(
        dict(energy=cse.energy, composition=formula), mp_elemental_ref_energies
    )
    for formula, cse in tqdm(
        df_eqnorm.set_index(Key.formula)[Key.computed_structure_entry].items(),
        total=len(df_eqnorm),
    )
]
df_eqnorm = df_eqnorm.round(4)

df_eqnorm.select_dtypes("number").to_csv(f"{out_path}.csv.gz")  # save csv storable
df_eqnorm.reset_index().to_json(
    f"{out_path}.json.gz", default_handler=as_dict_handler, orient="records", lines=True
)
