"""Module with code for calculating the thermal conductivity of a material.

Code is adapted from https://github.com/MPA2suite/k_SRME/blob/6ff4c867/k_srme/conductivity.py.
All credit to Balázs Póta, Paramvir Ahlawat, Gábor Csányi, Michele Simoncelli. See
https://arxiv.org/abs/2408.00755 for details.
It was ported to this repo in https://github.com/janosh/matbench-discovery/pull/196 to
implement parallelization across input structures which allows scaling thermal
conductivity metric to larger test sets.
"""

from typing import Any

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from fairchem.core.preprocessing.atoms_to_graphs import AtomsToGraphs
from phono3py.api_phono3py import Phono3py
from phonopy.structure.atoms import PhonopyAtoms
from torch_geometric.data import Batch
from tqdm import tqdm


def aseatoms2phonoatoms(atoms: Atoms) -> PhonopyAtoms:
    return PhonopyAtoms(
        atoms.symbols, cell=atoms.cell, positions=atoms.positions, pbc=True
    )


def aseatoms2phono3py(
    atoms: Atoms,
    fc2_supercell: np.typing.NDArray,
    fc3_supercell: np.typing.NDArray,
    primitive_matrix: np.typing.NDArray | None = None,
    **kwargs: Any,
) -> Phono3py:
    unitcell = aseatoms2phonoatoms(atoms)
    return Phono3py(
        unitcell=unitcell,
        supercell_matrix=fc3_supercell,
        phonon_supercell_matrix=fc2_supercell,
        primitive_matrix=primitive_matrix,
        **kwargs,
    )


def calculate_fc2_set(
    ph3: Phono3py, calculator: Calculator, pbar_kwargs: dict[str, Any] | None = None
) -> np.ndarray:
    """Calculate 2nd order force constants.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate force constants.
        calculator (Calculator): ASE calculator to compute forces.
        pbar_kwargs (dict[str, Any] | None): Arguments passed to tqdm progress bar.
            Defaults to None.

    Returns:
        np.ndarray: Array of forces for each displacement
    """

    forces = []
    n_atoms = len(ph3.phonon_supercell)

    displacements = ph3.phonon_supercells_with_displacements
    for supercell in tqdm(
        displacements,
        desc=f"FC2 calculation: {ph3.unitcell.formula}",
        **(pbar_kwargs or {}),
    ):
        if supercell is not None:
            atoms = Atoms(
                supercell.symbols,
                cell=supercell.cell,
                positions=supercell.positions,
                pbc=True,
            )
            atoms.calc = calculator
            force = atoms.get_forces()
        else:
            force = np.zeros((n_atoms, 3))
        forces += [force]

    force_set = np.array(forces)
    ph3.phonon_forces = force_set
    return force_set


def calculate_fc3_set(
    ph3: Phono3py,
    calculator: Calculator,
    pbar_kwargs: dict[str, Any] | None = None,
) -> np.ndarray:
    """Calculate 3rd order force constants.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate force constants.
        calculator (Calculator): ASE calculator to compute forces.
        pbar_kwargs (dict[str, Any] | None): Passed to tqdm progress bar.
            Defaults to None.

    Returns:
        np.ndarray: Array of forces for each displacement
    """
    forces: list[np.ndarray] = []
    n_atoms = len(ph3.supercell)

    desc = f"FC3 calculation: {ph3.unitcell.formula}"
    task_idx = (pbar_kwargs or {}).get("position")
    if task_idx:
        desc = f"{task_idx}. {desc}"
    displacements = ph3.supercells_with_displacements
    for supercell in tqdm(displacements, desc=desc, **pbar_kwargs or {}):
        if supercell is None:
            forces += [np.zeros((n_atoms, 3))]
        else:
            atoms = Atoms(
                supercell.symbols,
                cell=supercell.cell,
                positions=supercell.positions,
                pbc=True,
            )
            atoms.calc = calculator
            forces += [atoms.get_forces()]

    force_set = np.array(forces)
    ph3.forces = force_set
    return force_set


def init_phono3py(
    atoms: Atoms,
    *,
    log: bool = True,
    symprec: float = 1e-5,
    displacement_distance: float = 0.03,
    **kwargs: Any,
) -> tuple[Phono3py, list[Any], list[Any]]:
    """Calculate fc2 and fc3 force lists from phono3py.

    Args:


    Raises:


    Returns:

    """
    if not log:
        log_level = 0
    elif log is not None:
        log_level = 1

    formula = atoms.get_chemical_formula(mode="metal")
    if "fc2_supercell" not in atoms.info:
        raise ValueError(
            f'{formula} "fc2_supercell" was not found in atoms.info'
            "when calculating force sets."
        )

    if "fc3_supercell" not in atoms.info:
        raise ValueError(
            f'{formula} "fc3_supercell" was not found in atoms.info'
            "when calculating force sets."
        )

    if "fc3_supercell" not in atoms.info:
        raise ValueError(
            f'{formula} "q_mesh" was not found in atoms.info '
            "when calculating force sets."
        )

    # Initialise Phono3py object
    ph3 = aseatoms2phono3py(
        atoms,
        fc2_supercell=atoms.info["fc2_supercell"],
        fc3_supercell=atoms.info["fc3_supercell"],
        primitive_matrix="auto",
        symprec=symprec,
        log_level=log_level,
        **kwargs,
    )

    ph3.mesh_numbers = atoms.info["q_mesh"]

    ph3.generate_displacements(distance=displacement_distance)

    return ph3


def get_fc2_and_freqs(
    ph3: Phono3py, calculator: Calculator, pbar_kwargs: dict[str, Any] | None = None
) -> tuple[Phono3py, np.ndarray, np.ndarray]:
    """Calculate 2nd order force constants and phonon frequencies.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate force constants.
        calculator (Calculator): ASE calculator to compute forces.
        pbar_kwargs (dict[str, Any] | None): Arguments passed to tqdm progress bar.
            Defaults to None.

    Returns:
        tuple[Phono3py, np.ndarray, np.ndarray]: Tuple of (Phono3py object, force
            constants array, frequencies array [shape: (n_bz_grid, n_bands)])

    Raises:
        ValueError: If mesh_numbers not set
    """
    if ph3.mesh_numbers is None:
        raise ValueError(
            "mesh_numbers was not found in phono3py object and was not provided as "
            "an argument when calculating phonons from phono3py object."
        )

    pbar_kwargs = {"leave": False} | (pbar_kwargs or {})
    fc2_set = calculate_fc2_set(ph3, calculator, pbar_kwargs=pbar_kwargs)

    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.init_phph_interaction(symmetrize_fc3q=False)
    ph3.run_phonon_solver()

    freqs, _eigvecs, _grid = ph3.get_phonon_data()

    return ph3, fc2_set, freqs


def load_force_sets(
    ph3: Phono3py, fc2_set: np.ndarray, fc3_set: np.ndarray
) -> Phono3py:
    """Load pre-computed force sets into Phono3py object.

    Args:
        ph3 (Phono3py): Phono3py object to load force sets into
        fc2_set (np.ndarray): 2nd order force constants array
        fc3_set (np.ndarray): 3rd order force constants array

    Returns:
        Phono3py: Phono3py object with loaded force sets
    """
    ph3.phonon_forces = fc2_set
    ph3.forces = fc3_set
    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.produce_fc3(symmetrize_fc3r=True)

    return ph3


def calc_mode_kappa_tot(
    mode_kappa_p_rta: np.ndarray,
    mode_kappa_coherence: np.ndarray,
    heat_capacity: np.ndarray,
) -> np.ndarray:
    """Calculate total mode kappa from particle-like RTA and coherence terms.

    Args:
        mode_kappa_p_rta (np.ndarray): Mode kappa from particle-like RTA with shape
            (T, q-points, bands, xyz)
        mode_kappa_coherence (np.ndarray): Mode kappa from wave-like coherence with
            shape (T, q-points, bands, bands, xyz)
        heat_capacity (np.ndarray): Mode heat capacities with shape
            (T, q-points, bands)

    Returns:
        np.ndarray: Total (particle-like + wave-like) thermal conductivity per phonon
            mode with shape (T, q-points, bands)
    """
    # Temporarily silence divide warnings since we handle NaN values below
    with np.errstate(divide="ignore", invalid="ignore"):
        mode_kappa_c_per_mode = 2 * (  # None equiv to np.newaxis
            (mode_kappa_coherence * heat_capacity[:, :, :, None, None])
            / (heat_capacity[:, :, :, None, None] + heat_capacity[:, :, None, :, None])
        ).sum(axis=2)

    mode_kappa_c_per_mode[np.isnan(mode_kappa_c_per_mode)] = 0

    return mode_kappa_c_per_mode + mode_kappa_p_rta


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
        print("msg from robert.cho : Bug may happen do not trust the result ")
        g_idx = 0
        m_idx = 0
        while m_idx < len(multiply_0):
            if multiply_0[m_idx]:
                graph_list.insert(g_idx, graph_list[0].clone())
            m_idx = m_idx + 1
            g_idx = g_idx + 1
    forces = []
    import os

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
