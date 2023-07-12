# %%
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from pqdm.processes import pqdm
from pymatgen.core import Structure
from pymatgen.io.jarvis import JarvisAtomsAdaptor

from matbench_discovery import DEBUG
from matbench_discovery.data import DATA_FILES, df_wbm

__author__ = "Janosh Riebesell, Philipp Benner"
__date__ = "2023-06-03"


# %% read environment variables
batch = int(os.getenv("TASK_ID", default="0"))
out_dir = os.getenv("SBATCH_OUTPUT", default="2022-10-19-alignn-wbm-IS2RE")


# %%
n_splits = 100
n_processes_per_task = 10
module_dir = os.path.dirname(__file__)
# model_name = "mp_e_form_alignn"  # pre-trained by NIST
model_name = f"{module_dir}/data-train-result/best-model.pth"
task_type = "IS2RE"
target_col = "e_form_per_atom_mp2020_corrected"
input_col = "initial_structure"
id_col = "material_id"
job_name = f"{model_name}-wbm-{task_type}{'-debug' if DEBUG else ''}"
out_path = (
    f"{out_dir}/{'alignn-relaxed-structs' if batch == 0 else f'{batch=}'}.json.gz"
)

if batch < 0 or batch > n_splits:
    raise SystemExit(f"Invalid task_id={batch}")
if batch > 0 and not os.path.exists(out_dir):
    os.mkdir(out_dir)
if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exiting")


# %% Load data
data_path = {
    "IS2RE": DATA_FILES.wbm_initial_structures,
    "RS2RE": DATA_FILES.wbm_computed_structure_entries,
}[task_type]
input_col = {"IS2RE": "initial_structure", "RS2RE": "relaxed_structure"}[task_type]

df_in = pd.read_json(data_path).set_index(id_col)

df_in[target_col] = df_wbm[target_col]
if task_type == "RS2RE":
    df_in[input_col] = [x["structure"] for x in df_in.computed_structure_entry]
assert input_col in df_in, f"{input_col=} not in {list(df_in)}"

# Split data into parts and process only one batch
if batch != 0:
    df_in = np.array_split(df_in, 100)[batch - 1]
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
        device=f"cuda:{batch % 4}" if torch.cuda.is_available() else "cpu",
        logfile="/dev/null",
    )
    # Relax structure
    opt, _, _ = ff.optimize_atoms(trajectory=None, logfile="/dev/null")

    return JarvisAtomsAdaptor.get_structure(opt)


structures = [
    df_in.loc[material_id]["initial_structure"] for material_id in df_in.index
]
df_relaxed = pqdm(structures, alignn_relax, n_jobs=n_processes_per_task)

df_in = df_in.assign(relaxed_structure=df_relaxed)


# %% save results
df_in.to_json(out_path)

# Examples of materials that take ages to converge:
# task_id = 75, df_in.iloc[856]: wbm-3-76848
# task_id = 75, df_in.iloc[986]: wbm-3-76978
