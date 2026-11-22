"""Build the MPtrj per-frame summary cache read by scripts/export_data_fig_payloads.py
and paper/scripts/eda_mp_trj.py.
"""

# %%
from collections import defaultdict

import numpy as np
import pandas as pd
from kemio.io import Frame
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import MP_DIR
from matbench_discovery.data import frames_from_zip
from matbench_discovery.energy import formation_energies, mp_elemental_ref_energies
from matbench_discovery.enums import DataFiles, MbdKey

__author__ = "Janosh Riebesell"
__date__ = "2023-11-22"


mp_trj_summary_path = f"{MP_DIR}/2022-09-16-mp-trj-summary.json.bz2"


# %% read extXYZ frames from the zipped directory without unpacking the whole archive
# takes ~8 mins on M2 Max
# takes ~5 mins on M3 Max
frame_list = frames_from_zip(DataFiles.mp_trj_extxyz.path)

mp_trj_frames: dict[str, list[Frame]] = defaultdict(list)
for frame in frame_list:
    mp_id = frame.structure.info[Key.mat_id]
    assert mp_id.startswith(("mp-", "mvc-"))
    mp_trj_frames[mp_id].append(frame)

del frame_list  # free up memory

assert len(mp_trj_frames) == 145_923  # number of unique MP IDs


# %%
def info_dict_to_id(info: dict[str, int | str]) -> str:
    """Construct a unique frame ID from the frame info dict."""
    return f"{info[Key.task_id]}-{info['calc_id']}-{info['ionic_step']}"


df_mp_trj = pd.DataFrame(
    {
        info_dict_to_id(frame.structure.info): frame.structure.info
        | {
            Key.forces: frame.forces,
            # 3x3 Cauchy stress (eV/A^3), not the Voigt 6-vector of earlier caches
            Key.stress: frame.stress,
            # (n_atoms,) collinear moments like ASE's get_magnetic_moments()
            Key.magmoms: frame.structure.site_data.get("magmoms"),
            Key.formula: frame.structure.composition.formula,
            Key.atom_nums: frame.structure.symbols,
        }
        for frames in tqdm(mp_trj_frames.values(), total=len(mp_trj_frames))
        for frame in frames
    }
).T.convert_dtypes()  # convert object columns to float/int where possible
df_mp_trj.index.name = "frame_id"
assert len(df_mp_trj) == 1_580_395  # number of total frames
if Key.formula not in df_mp_trj:
    raise KeyError(f"{Key.formula!s} not in {df_mp_trj.columns=}")

# this is the unrelaxed (but MP2020 corrected) formation energy per atom of the actual
# relaxation step
df_mp_trj[MbdKey.e_form_dft] = formation_energies(
    df_mp_trj[Key.formula],
    df_mp_trj["mp2020_corrected_energy"].astype(float),
    mp_elemental_ref_energies,
)
df_mp_trj[Key.stress_trace] = [
    np.trace(stress) / 3 for stress in tqdm(df_mp_trj[Key.stress])
]


# %%
df_mp_trj.to_json(mp_trj_summary_path)
