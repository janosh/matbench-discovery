"""Build the MPtrj per-frame summary cache read by scripts/export_data_fig_payloads.py
and paper/scripts/eda_mp_trj.py.
"""

# %%
from collections import defaultdict

import ase
import ase.io.extxyz
import numpy as np
import pandas as pd
from pymatgen.core.tensors import Tensor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import MP_DIR
from matbench_discovery.data import ase_atoms_from_zip
from matbench_discovery.energy import calc_energy_from_e_refs, mp_elemental_ref_energies
from matbench_discovery.enums import DataFiles, MbdKey

__author__ = "Janosh Riebesell"
__date__ = "2023-11-22"


mp_trj_summary_path = f"{MP_DIR}/2022-09-16-mp-trj-summary.json.bz2"


# %% extract extXYZ files from zipped directory without unpacking the whole archive
# takes ~8 mins on M2 Max
# takes ~5 mins on M3 Max
atoms_list = ase_atoms_from_zip(DataFiles.mp_trj_extxyz.path)

mp_trj_atoms: dict[str, list[ase.Atoms]] = defaultdict(list)
for atoms in atoms_list:
    mp_id = atoms.info[Key.mat_id]
    assert mp_id.startswith(("mp-", "mvc-"))
    mp_trj_atoms[mp_id].append(atoms)

del atoms_list  # free up memory

assert len(mp_trj_atoms) == 145_923  # number of unique MP IDs


# %%
def info_dict_to_id(info: dict[str, int | str]) -> str:
    """Construct a unique frame ID from the atoms info dict."""
    return f"{info[Key.task_id]}-{info['calc_id']}-{info['ionic_step']}"


df_mp_trj = pd.DataFrame(
    {
        info_dict_to_id(atoms.info): atoms.info
        | {
            Key.forces: atoms.get_forces(),
            Key.stress: atoms.get_stress(),
            Key.magmoms: atoms.get_magnetic_moments()
            if "magmoms" in atoms.calc.results
            else None,
            Key.formula: str(atoms.symbols),
            Key.atom_nums: atoms.symbols,
        }
        for atoms_list in tqdm(mp_trj_atoms.values(), total=len(mp_trj_atoms))
        for atoms in atoms_list
    }
).T.convert_dtypes()  # convert object columns to float/int where possible
df_mp_trj.index.name = "frame_id"
assert len(df_mp_trj) == 1_580_395  # number of total frames
if Key.formula not in df_mp_trj:
    raise KeyError(f"{Key.formula!s} not in {df_mp_trj.columns=}")

# this is the unrelaxed (but MP2020 corrected) formation energy per atom of the actual
# relaxation step
df_mp_trj[MbdKey.e_form_dft] = [
    calc_energy_from_e_refs(
        {"composition": row[Key.formula], "energy": row["mp2020_corrected_energy"]},
        mp_elemental_ref_energies,
    )
    for _idx, row in tqdm(
        df_mp_trj.iterrows(), total=len(df_mp_trj), desc="Compute formation energies"
    )
]
df_mp_trj[Key.stress_trace] = [
    np.trace(Tensor.from_voigt(stress)) / 3 for stress in tqdm(df_mp_trj[Key.stress])
]


# %%
df_mp_trj.to_json(mp_trj_summary_path)
