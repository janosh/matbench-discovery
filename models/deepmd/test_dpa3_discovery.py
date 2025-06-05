"""Theoretically, one can use the test code of SevenNet or MACE to test a DP model.
For convenience, we used dflow to orchestrate the tests.
Below are the core functions.
"""

from __future__ import annotations

import pickle
from glob import glob
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE
from deepmd.calculator import DP
from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import as_dict_handler
from matbench_discovery.enums import Task

if TYPE_CHECKING:
    import numpy as np
    from ase import Atoms


class Relaxer:
    """Wrapper for ase.Atoms

    Args:
        model (Path): DP model, a path to the freezed model is needed.
    """

    def __init__(self, model: str | Path) -> None:
        try:
            self.calculator = DP(Path(model))
        except Exception as exc:
            print(f"DP calculator load failed: {exc}")

        self.optimizer = FIRE
        self.ase_adaptor = AseAtomsAdaptor()

    def relax(
        self, atoms: Atoms, fmax: float, steps: int, traj_file: str | None = None
    ) -> dict[str, Any]:
        """Relax atomic structure using ASE optimizer.

        Args:
            atoms (Atoms): Atomic structure to relax.
            fmax (float): Maximum force criterion for convergence.
            steps (int): Maximum number of optimization steps.
            traj_file (str | None, optional): Path to save trajectory. Defaults to None.

        Returns:
            dict[str, Any]: Dictionary containing final structure and trajectory.
        """
        if isinstance(atoms, Structure | Molecule):
            atoms = self.ase_adaptor.get_atoms(atoms)

        atoms.calc = self.calculator
        obs = TrajectoryObserver(atoms)
        atoms = FrechetCellFilter(atoms)
        opt = self.optimizer(atoms)
        opt.attach(obs)
        opt.run(fmax=fmax, steps=steps)
        obs()
        if traj_file is not None:
            obs.save(traj_file)
        if isinstance(atoms, FrechetCellFilter):
            atoms = atoms.atoms
        return {
            "final_structure": self.ase_adaptor.get_structure(atoms).as_dict(),
            "trajectory": obs,
        }


class TrajectoryObserver:
    """Trajectory observer is a hook in the relaxation process that saves the
    intermediate structures
    """

    def __init__(self, atoms: Atoms) -> None:
        """Initialize trajectory observer.

        Args:
            atoms (Atoms): The structure to observe.
        """
        self.atoms = atoms
        self.energies: list[float] = []
        self.forces: list[np.ndarray] = []
        self.stresses: list[np.ndarray] = []
        self.atom_positions: list[np.ndarray] = []
        self.cells: list[np.ndarray] = []

    def __call__(self) -> None:
        """Save properties of Atoms during relaxation."""
        self.energies.append(self.compute_energy())
        self.forces.append(self.atoms.get_forces())
        self.stresses.append(self.atoms.get_stress())
        self.atom_positions.append(self.atoms.get_positions())
        self.cells.append(self.atoms.get_cell()[:])

    def compute_energy(self) -> float:
        """Calculate the potential energy.

        Returns:
            float: Potential energy of the system.
        """
        return self.atoms.get_potential_energy()

    def save(self, filename: str) -> None:
        """Save the trajectory to file.

        Args:
            filename (str): Filename to save the trajectory.
        """
        with open(filename, mode="wb") as f:
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


def relax_run(
    filepath: str, model: str, relaxer: Relaxer, fmax: float = 0.05, steps: int = 500
) -> pd.DataFrame:
    """Run structure relaxation on input structures.

    Args:
        filepath (str): Path to input JSON file containing structures.
        model (str): Name of the model used for relaxation.
        relaxer (Relaxer): Relaxer instance to perform relaxations.
        fmax (float, optional): Force convergence criterion. Defaults to 0.05.
        steps (int, optional): Maximum optimization steps. Defaults to 500.

    Returns:
        pd.DataFrame: Results containing relaxed structures and energies.
    """
    task_type = Task.IS2RE

    df_in = pd.read_json(filepath)
    print("\nAll Data Loading Finished!!!\n")

    relax_results: dict[str, dict[str, Any]] = {}

    input_col = {Task.IS2RE: Key.initial_struct, Task.RS2RE: Key.final_struct}[
        task_type
    ]
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
    df_out.index.name = Key.mat_id
    return df_out


def relax_structures(
    input_dir: str, model: str | Path, out_path: str = "out.json.gz"
) -> None:
    """Relax structures from input directory using given model.

    Args:
        input_dir (str): Path to input JSON file.
        model (str | Path): Path to model file.
        out_path (str, optional): Path to output file. Defaults to "out.json.gz".
    """
    relaxer = Relaxer(model=model)

    ret_df = relax_run(input_dir, model="dp", relaxer=relaxer, fmax=0.05, steps=500)
    ret_df.reset_index().to_json(
        out_path, default_handler=as_dict_handler, orient="records", lines=True
    )


if __name__ == "__main__":
    input_dirs = sorted(glob("./data/wbm_data_*.json"))

    # this actually runs in parallel on multiple Nodes, orchestrated by dflow
    for input_dir in input_dirs:
        relax_structures(input_dir, "./dpa-3.1-3M-ft.pth")
