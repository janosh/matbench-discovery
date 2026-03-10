"""Batch-wise computation for third-order force constants (FC3)."""

import os
from typing import Any

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from fairchem.core.preprocessing import AtomsToGraphs
from phono3py.api_phono3py import Phono3py
from torch_geometric.data import Batch


def calculate_fc3_set_batch(
    ph3: Phono3py,
    calculator: Calculator,
    pbar_kwargs: dict[str, Any] | None = None,
) -> np.ndarray:
    """Calculate 3rd order force constants in batches for efficiency."""
    pbar_kwargs = pbar_kwargs or {}
    n_atoms = len(ph3.supercell)
    graph_list = []
    skip_supercell = []
    atoms_to_graphs = AtomsToGraphs(r_edges=False)

    for supercell in ph3.supercells_with_displacements:
        if supercell is not None:
            atoms = Atoms(
                supercell.symbols,
                cell=supercell.cell,
                positions=supercell.positions,
                pbc=True,
            )
            graph_list.append(atoms_to_graphs.convert(atoms))
            skip_supercell.append(False)
        else:
            skip_supercell.append(True)

    if len(skip_supercell) != len(graph_list):
        for idx, should_skip in enumerate(skip_supercell):
            if should_skip:
                graph_list.insert(idx, graph_list[0].clone())

    max_atoms_per_batch = int(os.environ.get("N_ATOMS", "128"))
    batch_size = max(max_atoms_per_batch // n_atoms, 1)
    device = "cuda"

    forces = []
    start_idx = 0
    while start_idx < len(graph_list):
        end_idx = min(start_idx + batch_size, len(graph_list))
        batch = Batch.from_data_list(graph_list[start_idx:end_idx]).to(device)
        result = calculator.trainer.predict(batch, per_image=False, disable_tqdm=True)
        forces.append(result["forces"].detach().cpu().numpy().reshape(-1, n_atoms, 3))
        start_idx = end_idx

    batch_force = np.concatenate(forces)
    batch_force[np.array(skip_supercell)] = 0
    return batch_force


def get_fc3_batch(
    ph3: Phono3py,
    calculator: Calculator,
    *,
    pbar_kwargs: dict[str, Any] | None = None,
) -> tuple[Phono3py, np.ndarray]:
    """Calculate 3rd order force constants using batched computation.

    Args:
        ph3: Phono3py object for which to calculate force constants.
        calculator: ASE calculator to compute forces.
        pbar_kwargs: Arguments passed to progress bar. Defaults to None.

    Returns:
        tuple[Phono3py, np.ndarray]: (Phono3py object, FC3 force set array)
    """
    pbar_kwargs = {"leave": False} | (pbar_kwargs or {})
    fc3_set = calculate_fc3_set_batch(ph3, calculator, pbar_kwargs=pbar_kwargs)
    return ph3, fc3_set
