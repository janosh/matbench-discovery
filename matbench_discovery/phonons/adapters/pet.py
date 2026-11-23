"""PET-specific symmetrized and batched force evaluation."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from matbench_discovery.phonons.adapters.standard import StandardKappaAdapter
from matbench_discovery.phonons.thermal_conductivity import (
    ase_batch_forces,
    batched_displacement_forces,
    fc2_frequencies,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from phono3py.api_phono3py import Phono3py

    from matbench_discovery.phonons.pipeline import KappaSettings


def pet_batch_forces(calculator: Calculator) -> Callable[[list[Atoms]], np.ndarray]:
    """Return a batch evaluator using the calculator's batched ``compute_energy`` API
    if it has one, else one ASE call per structure (current metatomic
    SymmetrizedCalculator versions, whose calls still batch their internal rotational
    quadrature).
    """
    compute_energy = getattr(calculator, "compute_energy", None)
    if not callable(compute_energy):
        return ase_batch_forces(calculator)
    return lambda batch_atoms: np.asarray(
        compute_energy(batch_atoms, compute_forces_and_stresses=True)["forces"]
    )


class PetKappaAdapter(StandardKappaAdapter):
    """Wrap PET for symmetry and batch its FC2/FC3 force evaluations."""

    name = "pet"

    def prepare_calculator(
        self, calculator: Calculator, settings: KappaSettings
    ) -> Calculator:
        """Wrap a metatomic calculator in its rotational symmetrizer."""
        from metatomic_ase import SymmetrizedCalculator

        if isinstance(calculator, SymmetrizedCalculator):
            return calculator
        # SymmetrizedCalculator subclasses ase Calculator, so no cast is needed
        return SymmetrizedCalculator(
            calculator, batch_size=settings.batch_size, include_inversion=False
        )

    def calculate_fc2(
        self,
        phono3py: Phono3py,
        calculator: Calculator,
        settings: KappaSettings,
    ) -> tuple[Phono3py, np.ndarray, np.ndarray]:
        """Calculate PET FC2 and frequencies with batched force calls."""
        force_set = batched_displacement_forces(
            phono3py.phonon_supercells_with_displacements,
            pet_batch_forces(calculator),
            batch_size=settings.batch_size,
            n_atoms=len(phono3py.phonon_supercell),
        )
        return phono3py, force_set, fc2_frequencies(phono3py, force_set)

    def calculate_fc3(
        self,
        phono3py: Phono3py,
        calculator: Calculator,
        settings: KappaSettings,
        *,
        max_evaluations: int | None = None,
    ) -> np.ndarray:
        """Calculate PET FC3 with batched force calls."""
        force_set = batched_displacement_forces(
            phono3py.supercells_with_displacements,
            pet_batch_forces(calculator),
            batch_size=settings.batch_size,
            n_atoms=len(phono3py.supercell),
            max_evaluations=max_evaluations,
        )
        phono3py.forces = force_set
        return force_set
