"""
Batchwise computation for FC3 computation

"""

import os
from typing import Any

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from fairchem.core.preprocessing.atoms_to_graphs import AtomsToGraphs
from phono3py.api_phono3py import Phono3py
from torch_geometric.data import Batch


def calculate_fc3_set_batch(
    ph3: Phono3py,
    calculator: Calculator,
    pbar_kwargs: dict[str, Any] | None = None,
) -> np.ndarray:
    # calculate FC3 force set
    if pbar_kwargs is None:
        pbar_kwargs = {}
    forces = []
    nat = len(ph3.supercell)
    graph_list = []
    multiply_0 = []
    a2g = AtomsToGraphs(r_edges=False)
    for sc in ph3.supercells_with_displacements:
        if sc is not None:
            atoms = Atoms(sc.symbols, cell=sc.cell, positions=sc.positions, pbc=True)

            graph = a2g.convert(atoms)
            graph_list.append(graph)
            multiply_0.append(False)
        else:
            multiply_0.append(True)
    if len(multiply_0) != len(graph_list):
        g_idx = 0
        m_idx = 0
        while m_idx < len(multiply_0):
            if multiply_0[m_idx]:
                graph_list.insert(g_idx, graph_list[0].clone())
            m_idx = m_idx + 1
            g_idx = g_idx + 1
    forces = []

    maximum_natom = int(os.environ.get("NATOMS", "128"))
    batchsize = max(maximum_natom // nat, 1)

    device = "cuda"

    sidx = 0
    while sidx < len(graph_list):
        eidx = min(sidx + batchsize, len(graph_list))
        batch = Batch.from_data_list(graph_list[sidx:eidx]).to(device)
        res = calculator.trainer.predict(batch, per_image=False, disable_tqdm=True)
        forces.append(res["forces"].detach().cpu().numpy().reshape(-1, nat, 3))
        sidx = eidx

    forces2 = forces

    batch_force = np.concatenate(forces2)
    idx = np.array(multiply_0)
    batch_force[idx, :, :] = 0
    return batch_force


def get_fc3_batch(
    ph3: Phono3py,
    calculator: Calculator,
    *,
    pbar_kwargs: dict[str, Any] | None = None,
) -> tuple[Phono3py, np.ndarray]:
    if pbar_kwargs is None:
        pbar_kwargs = {"leave": False}
    fc3_set = calculate_fc3_set_batch(ph3, calculator, pbar_kwargs=pbar_kwargs)
    return ph3, fc3_set
