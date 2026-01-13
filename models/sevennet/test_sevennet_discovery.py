"""Test SevenNet relaxation on the WBM dataset."""

# %%
import os
from typing import Any

import numpy as np
import pandas as pd
import torch
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from ase.optimize.optimize import Optimizer
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from sevenn.calculator import SevenNetCalculator
from tqdm import tqdm

from matbench_discovery import WBM_DIR, timestamp
from matbench_discovery.data import as_dict_handler, ase_atoms_from_zip
from matbench_discovery.enums import DataFiles, Task

__author__ = "Jinmu You"
__date__ = "2026-01-12"

"""History
2024-06-25 (Yutack Park): SevenNet-0, first version
2024-12-10 (Yutack Park): SevenNet-l3i5
2025-03-14 (Yutack Park): SevenNet-mf-ompa. Different variants can be selected.
2026-01-12 (Jinmu You): SevenNet-omni-i12 added.
"""


model_name = "sevennet"
model_variant = "sevennet-mf-ompa"  # choose 7net model variant to eval

device = "cuda" if torch.cuda.is_available() else "cpu"
calc_kwargs = {
    "sevennet-0": {"model": "7net-0"},
    "sevennet-l3i5": {"model": "7net-l3i5"},
    "sevennet-mf-ompa": {"model": "7net-mf-ompa", "modal": "mpa"},
    "sevennet-omni-i12": {"model": "7net-omni-i12", "modal": "mpa"},
}[model_variant]
calc_kwargs["device"] = device

# Will be removed after integrating model checkpoint download into sevenn package
checkpoint_urls = {
    "sevennet-omni-i12": "https://figshare.com/ndownloader/files/60977863",
}
if model_variant in checkpoint_urls:
    cache_dir = Path.home() / ".cache" / "sevennet"
    cache_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = cache_dir / f"checkpoint_{model_variant.replace('-', '_')}.pth"

    if not checkpoint_path.exists():
        print(f"Downloading {model_variant} checkpoint to {checkpoint_path}...")
        import requests
        response = requests.get(checkpoint_urls[model_variant], stream=True)
        response.raise_for_status()
        with open(checkpoint_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete.")
    else:
        print(f"Using cached checkpoint: {checkpoint_path}")
    calc_kwargs["model"] = str(checkpoint_path)

# %% this config is editable
smoke_test = True
model_name = model_variant
task_type = Task.IS2RE
job_name = f"{model_name}-wbm-{task_type}"
ase_optimizer = "FIRE"
device = "cuda" if torch.cuda.is_available() else "cpu"

# They gives almost the same results. These values are for reproducibility.
max_steps = 800 if model_variant in ["sevennet-mf-ompa", "sevennet-omni-i12"] else 500
force_max = 0.02 if model_variant in ["sevennet-mf-ompa", "sevennet-omni-i12"] else 0.05

slurm_array_task_count = 32


# %%
slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", "debug")

os.makedirs(out_dir := "./results", exist_ok=True)
out_path = f"{out_dir}/{model_name}-{slurm_array_task_id:>03}.json.gz"

data_path = {Task.IS2RE: DataFiles.wbm_initial_atoms.path}[task_type]
print(f"\nJob {job_name!r} running {timestamp}", flush=True)
print(f"{data_path=}", flush=True)

# Initialize ASE SevenNet Calculator from checkpoint
seven_net_calc = SevenNetCalculator(**calc_kwargs)


# %%
print(f"Read data from {data_path}")
zip_filename = f"{WBM_DIR}/2024-08-04-wbm-initial-atoms.extxyz.zip"
atoms_list = ase_atoms_from_zip(zip_filename)

if slurm_array_job_id == "debug":
    if smoke_test:
        atoms_list = atoms_list[:128]
    else:
        pass
elif slurm_array_task_count > 1:
    atoms_list = np.array_split(atoms_list, slurm_array_task_count)[
        slurm_array_task_id - 1
    ]

relax_results: dict[str, dict[str, Any]] = {}
optim_cls: type[Optimizer] = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]


# %%
for atoms in tqdm(atoms_list, desc="Relaxing"):
    mat_id = atoms.info[Key.mat_id]
    if mat_id in relax_results:
        continue
    try:
        atoms.calc = seven_net_calc
        if max_steps > 0:
            atoms = FrechetCellFilter(atoms)
            optimizer = optim_cls(atoms, logfile=None)
            optimizer.run(fmax=force_max, steps=max_steps)
        energy = atoms.get_potential_energy()  # relaxed energy
        # if max_steps > 0, atoms is wrapped by FrechetCellFilter, so need to getattr
        unwrapped = atoms.atoms if hasattr(atoms, "atoms") else atoms
        relaxed_struct = AseAtomsAdaptor.get_structure(unwrapped)
        relax_results[mat_id] = {"structure": relaxed_struct, "energy": energy}
    except Exception as exc:
        print(f"Failed to relax {mat_id}: {exc!r}")
        continue

df_out = pd.DataFrame(relax_results).T.add_prefix("sevennet_")
df_out.index.name = Key.mat_id


# %%
if not smoke_test:
    df_out.reset_index().to_json(
        out_path, default_handler=as_dict_handler, orient="records", lines=True
    )