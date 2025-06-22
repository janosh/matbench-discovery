"""Test Eqnorm relaxation on the WBM dataset.

History:
- 2025-03-28 (Yuzhuo Chen): Eqnorm Mptrj, first version
"""

# %%
import os
from typing import Any

import pandas as pd
import torch
from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from eqnorm.calculator import EqnormCalculator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import WBM_DIR, timestamp
from matbench_discovery.data import as_dict_handler, ase_atoms_from_zip
from matbench_discovery.enums import DataFiles, Task


def process_and_save(atoms_list: list[Atoms], out_dir: str, job_id: int) -> None:
    optim_cls = {"FIRE": FIRE, "LBFGS": LBFGS}[ase_optimizer]
    out_path = f"{out_dir}/{model_name}-{job_id:>03}.json.gz"

    relax_results: dict[str, dict[str, Any]] = {}
    for atoms in tqdm(atoms_list, desc="Relaxing"):
        mat_id = atoms.info[Key.mat_id]
        if mat_id in relax_results:
            continue
        try:
            atoms.calc = calc
            if max_steps > 0:
                atoms = FrechetCellFilter(atoms)
                optimizer = optim_cls(atoms, logfile="/dev/null")
                optimizer.run(fmax=force_max, steps=max_steps)
            energy = atoms.get_potential_energy()  # relaxed energy
            # if max_steps > 0, atoms is wrapped by FrechetCellFilter
            relaxed_struct = AseAtomsAdaptor.get_structure(
                getattr(atoms, "atoms", atoms)
            )
            relax_results[mat_id] = {"structure": relaxed_struct, "energy": energy}
        except Exception as exc:
            print(f"Failed to relax {mat_id}: {exc!r}")
            continue

    df_out = pd.DataFrame(relax_results).T.add_prefix(f"{model_name}_")
    df_out.index.name = Key.mat_id

    # %%
    df_out.reset_index().to_json(
        out_path, default_handler=as_dict_handler, orient="records", lines=True
    )


if __name__ == "__main__":
    model_name = "eqnorm"
    model_variant = "eqnorm-mptrj"
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # device = torch.device("cpu")
    ase_optimizer = "FIRE"
    out_dir = "./results"
    os.makedirs(out_dir, exist_ok=True)

    max_steps = 500
    force_max = 0.02

    slurm_array_task_count = int(os.getenv("TOTAL_TASKS", "1"))
    slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))

    task_type = Task.IS2RE
    job_name = f"{model_name}-wbm-{task_type}"
    data_path = {Task.IS2RE: DataFiles.wbm_initial_atoms.path}[task_type]
    print(f"\nJob {job_name!r} running {timestamp}", flush=True)
    print(f"{data_path=}", flush=True)

    print(f"Read data from {data_path}")
    zip_filename = f"{WBM_DIR}/2024-08-04-wbm-initial-atoms.extxyz.zip"
    atoms_list = ase_atoms_from_zip(zip_filename)

    atoms_list = atoms_list[slurm_array_task_id::slurm_array_task_count]

    # load my own calculator
    calc = EqnormCalculator(
        model_name=model_name,
        model_variant=model_variant,
        device=device,
    )

    process_and_save(atoms_list, out_dir, slurm_array_task_id)
