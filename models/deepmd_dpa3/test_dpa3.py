# Theoritically, one can use the test code of sevennet or mace to test a DP model. For convinience, we used dflow to orchestrate the tests, below are the core functions.

from __future__ import annotations

from matbench_discovery.enums import Key, Task
import ase
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
from pymatgen.core.structure import Molecule, Structure
from deepmd.calculator import DP as DPCalculator
from ase.optimize import FIRE

import torch
from ase.constraints import ExpCellFilter
import pickle
from typing import Optional, Union
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from matbench_discovery.data import as_dict_handler

OPTIMIZERS = {
    "FIRE": FIRE,
}


class Relaxer:
    """Wrapper for ase.Atoms

    Parameters:
    ----------
    model: Union[Path, str]
        Indicates which calculator to use during relaxation, `mace` will call default MACE-medium model,
        for DP model, a path to the freezed model is needed.
    optimizer: str
        The optimizer from ASE, supports `FIRE`, 
    relax_cell: bool
        Whether to relax cell with `ExpCellFilter`.
    """

    def __init__(
        self,
        model: Union[str, Path],
        optimizer: Optional[str] = "FIRE",
        relax_cell: Optional[bool] = True,
    ):
        if isinstance(model, Path):
            try:
                self.calculator = DPCalculator(model)
            except Exception as e:
                raise ValueError(f"DP calculator load failed: {e}")
        else:
            raise NotImplementedError(
                "Only DP calculators are supported."
            )


        self.optimizer = OPTIMIZERS[optimizer]
        self.relax_cell = relax_cell
        self.ase_adaptor = AseAtomsAdaptor()

    def relax(self, atoms, fmax: float, steps: int, traj_file: str = None):
        if isinstance(atoms, (Structure, Molecule)):
            atoms = self.ase_adaptor.get_atoms(atoms)

        atoms.set_calculator(self.calculator)
        obs = TrajectoryObserver(atoms)
        if self.relax_cell:
            atoms = ExpCellFilter(atoms)
        opt = self.optimizer(atoms)
        opt.attach(obs)
        opt.run(fmax=fmax, steps=steps)
        obs()
        if traj_file is not None:
            obs.save(traj_file)
        if isinstance(atoms, ExpCellFilter):
            atoms = atoms.atoms
        return {
            "final_structure": self.ase_adaptor.get_structure(atoms).as_dict(),
            "trajectory": obs,
        }

class TrajectoryObserver:
    """
    Trajectory observer is a hook in the relaxation process that saves the
    intermediate structures
    """

    def __init__(self, atoms: ase.Atoms):
        """
        Args:
            atoms (Atoms): the structure to observe
        """
        self.atoms = atoms
        self.energies: list[float] = []
        self.forces: list[np.ndarray] = []
        self.stresses: list[np.ndarray] = []
        self.atom_positions: list[np.ndarray] = []
        self.cells: list[np.ndarray] = []

    def __call__(self):
        """
        The logic for saving the properties of an Atoms during the relaxation
        Returns:
        """
        self.energies.append(self.compute_energy())
        self.forces.append(self.atoms.get_forces())
        self.stresses.append(self.atoms.get_stress())
        self.atom_positions.append(self.atoms.get_positions())
        self.cells.append(self.atoms.get_cell()[:])

    def compute_energy(self) -> float:
        """
        calculate the energy, here we just use the potential energy
        Returns:
        """
        energy = self.atoms.get_potential_energy()
        return energy

    def save(self, filename: str):
        """
        Save the trajectory to file
        Args:
            filename (str): filename to save the trajectory
        Returns:
        """
        with open(filename, "wb") as f:
            pickle.dump(
                {
                    "energy": self.energies,
                    "forces": self.forces,
                    "stresses": self.stresses,
                    "atom_positions": self.atom_positions,
                    "cell": self.cells,
                    "atomic_number": self.atoms.get_atomic_numbers(),
                },
                f,
            )

def relax_run(fpth: str, model: str, relaxer: Relaxer, fmax: float = 0.05, steps: int = 500):
    task_type = Task.IS2RE

    df_in = pd.read_json(fpth)
    print("\nAll Data Loading Finished!!!\n")

    relax_results: dict[str, dict[str, Any]] = {}

    input_col = {Task.IS2RE: Key.init_struct, Task.RS2RE: Key.final_struct}[task_type]
    structures = df_in[input_col].map(Structure.from_dict).to_dict()

    for material_id in tqdm(structures, desc="Relaxing"):
        if material_id in relax_results:
            continue
        try:
            result = relaxer.relax(structures[material_id], fmax=fmax, steps=steps)
            relax_results[material_id] = {
                f"{model}_structure": result["final_structure"],
                f"{model}_energy": result["trajectory"].energies[-1],
            }

        except Exception as exc:
            print(f"Failed to relax {material_id}: {exc!r}")

    df_out = pd.DataFrame(relax_results).T
    print("\nSaved to df.\n")
    df_out.index.name = Key.mat_id
    return df_out

import os
from pathlib import Path
from dflow import Step, Workflow, upload_artifact
from dflow.plugins.dispatcher import DispatcherExecutor
from dflow.python import OP, Artifact, PythonOPTemplate, Slices


def relax_structures(inputs, model) -> {"res"}:
    
    relaxer = Relaxer(model=model, optimizer="FIRE")

    ret_df = relax_run(inputs, model="dp", relaxer=relaxer, fmax=0.05, steps=500)
    ret_df.reset_index().to_json("out.json.gz", default_handler=as_dict_handler)
    return {"res": Path("out.json.gz")}


if __name__ == "__main__":
    input_dirs = ["./data/wbm_data_0.json", "./data/wbm_data_1.json", "./data/wbm_data_2.json", "./data/wbm_data_3.json", "./data/wbm_data_4.json", "./data/wbm_data_5.json",
                "./data/wbm_data_6.json", "./data/wbm_data_7.json", "./data/wbm_data_8.json", "./data/wbm_data_9.json", "./data/wbm_data_10.json", "./data/wbm_data_11.json",
                "./data/wbm_data_12.json", "./data/wbm_data_13.json", "./data/wbm_data_14.json", "./data/wbm_data_15.json", "./data/wbm_data_16.json", "./data/wbm_data_17.json",
                "./data/wbm_data_18.json", "./data/wbm_data_19.json", "./data/wbm_data_20.json", "./data/wbm_data_21.json", "./data/wbm_data_22.json", "./data/wbm_data_23.json",
                "./data/wbm_data_24.json","./data/wbm_data_25.json","./data/wbm_data_26.json",
                 "./data/wbm_data_27.json","./data/wbm_data_28.json", "./data/wbm_data_29.json",
                 "./data/wbm_data_30.json","./data/wbm_data_31.json","./data/wbm_data_32.json",
                 "./data/wbm_data_33.json","./data/wbm_data_34.json","./data/wbm_data_35.json",
                 "./data/wbm_data_36.json","./data/wbm_data_37.json",
                 "./data/wbm_data_38.json","./data/wbm_data_39.json"]

    relax_structures(input_dirs, Path("./2025-01-10-dpa3-openlam.pth"))

