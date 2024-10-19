# %%
import os

import numpy as np
import pandas as pd
from pymatgen.core import Structure
from pymatgen.io.jarvis import JarvisAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import today
from matbench_discovery.data import DataFiles, df_wbm
from matbench_discovery.enums import MbdKey, Task

__author__ = "Janosh Riebesell, Philipp Benner"
__date__ = "2023-07-11"


# %% read environment variables
task_id = int(os.getenv("TASK_ID", default="0"))
out_dir = os.getenv("SBATCH_OUTPUT", default=f"{today}-alignn-wbm-IS2RE")


# %%
n_splits = 100
n_processes_per_task = 10
module_dir = os.path.dirname(__file__)
# model_name = "mp_e_form_alignn"  # pre-trained by NIST
model_name = f"{out_dir}/best-model.pth"
task_type = Task.IS2RE
input_col = Key.init_struct
job_name = f"{model_name}-wbm-{task_type}"
out_path = (
    f"{out_dir}/{'alignn-relaxed-structs' if task_id == 0 else f'{task_id=}'}.json.gz"
)

if not (0 <= task_id <= n_splits):
    raise SystemExit(f"Invalid {task_id=}")
if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exiting")
os.makedirs(out_dir, exist_ok=True)


# %% Load data
data_path = {
    Task.IS2RE: DataFiles.wbm_initial_structures.path,
    Task.RS2RE: DataFiles.wbm_computed_structure_entries.path,
}[task_type]
input_col = {Task.IS2RE: Key.init_struct, Task.RS2RE: Key.final_struct}[task_type]

df_in = pd.read_json(data_path).set_index(Key.mat_id)

df_in[MbdKey.e_form_dft] = df_wbm[MbdKey.e_form_dft]
if task_type == Task.RS2RE:
    df_in[input_col] = [cse["structure"] for cse in df_in[Key.cse]]
if input_col not in df_in:
    raise TypeError(f"{input_col} not in {df_in.columns=}")

# Split data into parts and process only one batch
if task_id != 0:
    df_in = np.array_split(df_in, 100)[task_id - 1]
    print(f"Relaxing materials in range {df_in.index[0]} - {df_in.index[-1]}")
else:
    print("Relaxing full range of materials")


# %% Relax structures
def alignn_relax(structure: Structure) -> Structure:
    """Relax structure using Alignn FF.

    Args:
        structure (Structure): pymatgen object to relax.

    Returns:
        Structure: Relaxed structure.
    """
    # Cuda must be only initialized in child processes
    import torch
    from alignn.ff.ff import ForceField, default_path

    ff = ForceField(
        jarvis_atoms=JarvisAtomsAdaptor.get_atoms(Structure.from_dict(structure)),
        model_path=default_path(),
        device=f"cuda:{task_id % 4}" if torch.cuda.is_available() else "cpu",
        logfile="/dev/null",
    )
    # Relax structure
    opt, _, _ = ff.optimize_atoms(trajectory=None, logfile="/dev/null")

    return JarvisAtomsAdaptor.get_structure(opt)


structures = [df_in.loc[material_id][Key.init_struct] for material_id in df_in.index]
df_relaxed = tqdm(structures, alignn_relax, n_jobs=n_processes_per_task)

df_in = df_in.assign(relaxed_structure=df_relaxed)


# %% save results
df_in.to_json(out_path)

# Examples of materials that take ages to converge:
# task_id = 75, df_in.iloc[856]: wbm-3-76848
# task_id = 75, df_in.iloc[986]: wbm-3-76978
