# /// script
# requires-python = ">=3.11,<3.13"
# dependencies = [
#     "tqdm",
#     "matbench-discovery",
#     "pandas",
#     "pymatgen",
#     "pymatviz",
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../../", editable = true }
# ///

# modified from eqnorm script

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

e_form_nequix_col = "e_form_per_atom_nequix"
results = "./results"
pot_name = "nequix"
out_path = f"{results}/{pot_name}"
files = sorted(glob(f"{results}/{pot_name}-*.json.gz"))

dfs = {}
for file_path in tqdm(files):
    if file_path in dfs:
        continue
    df_i = pd.read_json(file_path, lines=True).set_index(Key.mat_id)
    dfs[file_path] = df_i

df_nequix = pd.concat(dfs.values()).round(4)

if len(df_nequix) != len(df_wbm):  # make sure there is no missing structure
    warnings.warn("Missing structures in nequix results", stacklevel=2)

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
    df_nequix.itertuples(), total=len(df_nequix), desc="ML energies to CSEs"
):
    struct_dict = getattr(row, f"{pot_name}_structure")
    energy = getattr(row, f"{pot_name}_energy")
    mlip_struct = Structure.from_dict(struct_dict)
    cse = df_wbm_cse.loc[row.Index, Key.computed_structure_entry]
    cse._energy = energy  # cse._energy is the uncorrected energy from MPtrj dataset (or vasp raw)  # noqa: E501, SLF001
    cse._structure = mlip_struct  # noqa: SLF001
    df_nequix.loc[row.Index, Key.computed_structure_entry] = cse


orig_len = len(df_nequix)
MaterialsProject2020Compatibility().process_entries(  # change in-place
    df_nequix[Key.computed_structure_entry], verbose=True, clean=True
)
if len(df_nequix) != orig_len:
    warnings.warn("Some structures were removed during energy correction", stacklevel=2)

df_nequix[Key.formula] = df_wbm[Key.formula]

df_nequix[e_form_nequix_col] = [
    # see https://matbench-discovery.materialsproject.org/data#mp-elemental-reference-energies
    # MP ref energies are the lowest energies found for unary structures of each element
    get_e_form_per_atom(
        dict(energy=cse.energy, composition=formula), mp_elemental_ref_energies
    )
    for formula, cse in tqdm(
        df_nequix.set_index(Key.formula)[Key.computed_structure_entry].items(),
        total=len(df_nequix),
    )
]
df_nequix = df_nequix.round(4)

df_nequix.select_dtypes("number").to_csv(f"{out_path}.csv.gz")  # save csv storable
df_nequix.reset_index().to_json(
    f"{out_path}.jsonl.gz",
    default_handler=as_dict_handler,
    orient="records",
    lines=True,
)
