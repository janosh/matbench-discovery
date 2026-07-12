"""PET-specific symmetrized and batched force evaluation."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

import numpy as np

from matbench_discovery.phonons.adapters.standard import StandardKappaAdapter
from matbench_discovery.phonons.thermal_conductivity import batched_displacement_forces

if TYPE_CHECKING:
    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from phono3py.api_phono3py import Phono3py
    from phonopy.structure.atoms import PhonopyAtoms

    from matbench_discovery.phonons.pipeline import KappaSettings


def _batched_force_set(
    displacements: list[PhonopyAtoms | None],
    calculator: Calculator,
    *,
    batch_size: int,
    n_atoms: int,
    max_evaluations: int | None = None,
) -> np.ndarray:
    """Evaluate displaced supercells in batches while retaining null entries."""

    def evaluate_batch(batch_atoms: list[Atoms]) -> np.ndarray:
        """Evaluate one PET displacement batch through its available API."""
        compute_energy = getattr(calculator, "compute_energy", None)
        if callable(compute_energy):
            outputs = compute_energy(batch_atoms, compute_forces_and_stresses=True)
            return np.asarray(outputs["forces"])
        # Current metatomic SymmetrizedCalculator exposes only ASE's single-system
        # API; each call still batches its internal rotational quadrature.
        for displaced_atoms in batch_atoms:
            displaced_atoms.calc = calculator
        return np.asarray([atoms.get_forces() for atoms in batch_atoms])

    return batched_displacement_forces(
        displacements,
        evaluate_batch,
        batch_size=batch_size,
        n_atoms=n_atoms,
        max_evaluations=max_evaluations,
    )


class PetKappaAdapter(StandardKappaAdapter):
    """Wrap PET for symmetry and batch its FC2/FC3 force evaluations."""

    name = "pet"

    def prepare_calculator(
        self, calculator: Calculator, settings: KappaSettings
    ) -> Calculator:
        """Wrap a metatomic calculator in its rotational symmetrizer."""
        from metatomic.torch.ase_calculator import SymmetrizedCalculator

        if isinstance(calculator, SymmetrizedCalculator):
            return calculator
        return cast(
            "Calculator",
            SymmetrizedCalculator(
                calculator,
                batch_size=settings.batch_size,
                include_inversion=False,
            ),
        )

    def calculate_fc2(
        self,
        phono3py: Phono3py,
        calculator: Calculator,
        settings: KappaSettings,
        *,
        progress: dict[str, Any] | None = None,  # noqa: ARG002
    ) -> tuple[Phono3py, np.ndarray, np.ndarray]:
        """Calculate PET FC2 and frequencies with batched force calls."""
        force_set = _batched_force_set(
            list(phono3py.phonon_supercells_with_displacements),
            calculator,
            batch_size=settings.batch_size,
            n_atoms=len(phono3py.phonon_supercell),
        )
        phono3py.phonon_forces = force_set
        phono3py.produce_fc2(symmetrize_fc2=True)
        phono3py.init_phph_interaction(symmetrize_fc3q=False)
        phono3py.run_phonon_solver()
        frequencies, _eigenvectors, _grid = phono3py.get_phonon_data()
        return phono3py, force_set, frequencies

    def calculate_fc3(
        self,
        phono3py: Phono3py,
        calculator: Calculator,
        settings: KappaSettings,
        *,
        progress: dict[str, Any] | None = None,  # noqa: ARG002
        max_evaluations: int | None = None,
    ) -> np.ndarray:
        """Calculate PET FC3 with batched force calls."""
        force_set = _batched_force_set(
            list(phono3py.supercells_with_displacements),
            calculator,
            batch_size=settings.batch_size,
            n_atoms=len(phono3py.supercell),
            max_evaluations=max_evaluations,
        )
        phono3py.forces = force_set
        return force_set
