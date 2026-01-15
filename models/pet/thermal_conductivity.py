from typing import Any

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from phono3py.api_phono3py import Phono3py
from phonopy.structure.atoms import PhonopyAtoms
from tqdm.auto import tqdm


def init_phono3py(
    atoms: Atoms,
    *,
    fc2_supercell: np.ndarray,
    fc3_supercell: np.ndarray,
    q_point_mesh: tuple[int, int, int] = (20, 20, 20),
    displacement_distance: float = 0.01,
    symprec: float = 1e-5,
    is_plusminus: bool = False,
    **kwargs: Any,
) -> Phono3py:
    """Initialize Phono3py object from ASE Atoms.

    Args:
        atoms (Atoms): ASE Atoms object to initialize from.
        fc2_supercell (np.ndarray): Supercell matrix for 2nd order force constants.
        fc3_supercell (np.ndarray): Supercell matrix for 3rd order force constants.
        q_point_mesh (tuple[int, int, int]): Mesh size for q-point sampling. Defaults
            to (20, 20, 20).
        displacement_distance (float): Displacement distance for force calculations.
            Defaults to 0.01.
        symprec (float): Symmetry precision for finding space group. Defaults to 1e-5.
        **kwargs (Any): Passed to Phono3py constructor.

    Returns:
        Phono3py: Initialized Phono3py object

    Raises:
        ValueError: If required metadata is missing from atoms.info
    """
    unit_cell = PhonopyAtoms(atoms.symbols, cell=atoms.cell, positions=atoms.positions)
    ph3 = Phono3py(
        unitcell=unit_cell,
        supercell_matrix=fc3_supercell,
        phonon_supercell_matrix=fc2_supercell,
        primitive_matrix="auto",
        symprec=symprec,
        **kwargs,
    )
    ph3.mesh_numbers = q_point_mesh

    ph3.generate_displacements(
        distance=displacement_distance, is_plusminus=is_plusminus
    )

    return ph3


def calculate_fc2_set(
    ph3: Phono3py,
    calculator: Calculator,
    pbar_kwargs: dict[str, Any] | None = None,
    batch_size: int = 1,
) -> np.ndarray:
    """Calculate 2nd order force constants. Requires initializing Phono3py with an FC2
    supercell matrix.

    Args:
        ph3 (Phono3py): Phono3py object for which to calculate force constants.
        calculator (Calculator): ASE calculator to compute forces.
        pbar_kwargs (dict[str, Any] | None): Arguments passed to tqdm progress bar.
            Defaults to None.

    Returns:
        np.ndarray: Array of forces for each displacement
    """
    print(f"Computing FC2 force set in {ph3.unitcell.formula}.")

    forces: list[np.ndarray] = []

    displacements = ph3.phonon_supercells_with_displacements
    for batch in tqdm(
        range(0, len(displacements), batch_size),
        desc=f"FC2 calculation: {ph3.unitcell.formula}",
        **pbar_kwargs or {},
    ):
        batch_displacements = displacements[batch : batch + batch_size]
        batch_atoms = [
            Atoms(
                supercell.symbols,
                cell=supercell.cell,
                positions=supercell.positions,
                pbc=True,
            )
            for supercell in batch_displacements
        ]
        res = calculator.compute_energy(  # type: ignore[attr-defined]
            batch_atoms, compute_forces_and_stresses=True
        )
        f = res["forces"]
        forces.extend(f)

    force_set = np.stack(forces)
    ph3.phonon_forces = force_set
    return force_set


def calculate_fc3_set(
    ph3: Phono3py,
    calculator: Calculator,
    pbar_kwargs: dict[str, Any] | None = None,
    batch_size: int = 1,
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

    desc = f"FC3 calculation: {ph3.unitcell.formula}"
    task_idx = (pbar_kwargs or {}).get("position")
    if task_idx:
        desc = f"{task_idx}. {desc}"
    displacements = ph3.supercells_with_displacements
    for batch in tqdm(
        range(0, len(displacements), batch_size),
        desc=f"FC3 calculation: {ph3.unitcell.formula}",
        **pbar_kwargs or {},
    ):
        batch_displacements = displacements[batch : batch + batch_size]
        batch_atoms = [
            Atoms(
                supercell.symbols,
                cell=supercell.cell,
                positions=supercell.positions,
                pbc=True,
            )
            for supercell in batch_displacements
        ]
        res = calculator.compute_energy(  # type: ignore[attr-defined]
            batch_atoms, compute_forces_and_stresses=True
        )
        f = res["forces"]
        forces.extend(f)

    force_set = np.stack(forces)
    ph3.forces = force_set
    return force_set


def get_fc2_and_freqs(
    ph3: Phono3py,
    calculator: Calculator,
    pbar_kwargs: dict[str, Any] | None = None,
    batch_size: int = 1,
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
    fc2_set = calculate_fc2_set(
        ph3,
        calculator,
        pbar_kwargs={"leave": False} | ({"disable": True}),
        batch_size=batch_size,
    )

    ph3.produce_fc2(symmetrize_fc2=True)
    ph3.init_phph_interaction(symmetrize_fc3q=False)
    ph3.run_phonon_solver()

    freqs, _eigvecs, _grid = ph3.get_phonon_data()

    return ph3, fc2_set, freqs
