# /// script
# requires-python = ">=3.11,<3.13"
# dependencies = [
#     "ase",
#     "matbench-discovery",
#     "nequix",
#     "pandas",
#     "pymatgen",
#     "pymatviz",
# ]
#
# [tool.uv.sources]
# nequix = { git = "https://github.com/atomicarchitects/nequix" }
# matbench-discovery = { path = "../../", editable = true }
# ///

# modified from eqnorm script

import os
from typing import Any

import pandas as pd
from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from nequix.calculator import NequixCalculator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import WBM_DIR
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

    df_out.reset_index().to_json(
        out_path, default_handler=as_dict_handler, orient="records", lines=True
    )


if __name__ == "__main__":
    model_name = "nequix"
    task_type = Task.IS2RE
    ase_optimizer = "FIRE"
    max_steps = 500
    force_max = 0.02
    out_dir = "./results"
    os.makedirs(out_dir, exist_ok=True)

    slurm_array_task_count = int(os.getenv("SLURM_ARRAY_TASK_COUNT", "1"))
    slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
    print(f"slurm_array_task_id: {slurm_array_task_id}")
    print(f"slurm_array_task_count: {slurm_array_task_count}")

    job_name = f"{model_name}-wbm-{task_type}"
    data_path = DataFiles.wbm_initial_atoms.path
    zip_filename = f"{WBM_DIR}/2024-08-04-wbm-initial-atoms.extxyz.zip"
    atoms_list = ase_atoms_from_zip(zip_filename)
    atoms_list = atoms_list[slurm_array_task_id::slurm_array_task_count]

    calc = NequixCalculator("nequix-mp-1")
    process_and_save(atoms_list, out_dir, slurm_array_task_id)
