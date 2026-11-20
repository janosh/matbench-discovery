"""PET-specific symmetrized and batched force evaluation."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np

from matbench_discovery.phonons.adapters.standard import StandardKappaAdapter

if TYPE_CHECKING:
    from matscikit import Structure
    from matscikit.model import Model as ForceModel

    from matbench_discovery.phonons.pipeline import KappaSettings


class PetBatchedModel:
    """PET behind the matscikit model protocol, batching through ``compute_energy``.

    Current metatomic SymmetrizedCalculator versions without ``compute_energy`` fall
    back to the single-structure ASE API; each call still batches its internal
    rotational quadrature.
    """

    def __init__(self, calculator: object) -> None:
        """Wrap a (symmetrized) metatomic ASE calculator."""
        self.calculator = calculator

    def compute(
        self, structures: Sequence[Structure], properties: Sequence[str]
    ) -> list[dict[str, Any]]:
        """Evaluate a batch of structures."""
        from matscikit.interop import AseCalculatorModel, structure_to_ase

        compute_energy = getattr(self.calculator, "compute_energy", None)
        if not callable(compute_energy) or set(properties) - {"energy", "forces"}:
            return AseCalculatorModel(self.calculator).compute(structures, properties)
        outputs = compute_energy(
            [structure_to_ase(structure) for structure in structures],
            compute_forces_and_stresses=True,
        )
        # forces stay per structure: batches may mix atom counts (ragged arrays)
        energies, forces = np.asarray(outputs["energy"]), outputs["forces"]
        return [
            {"energy": float(energies[idx]), "forces": np.asarray(forces[idx])}
            for idx in range(len(structures))
        ]


class PetKappaAdapter(StandardKappaAdapter):
    """Wrap PET for symmetry and batch its FC2/FC3 force evaluations."""

    name = "pet"

    def prepare_model(self, model: ForceModel, settings: KappaSettings) -> ForceModel:
        """Wrap the metatomic calculator in its rotational symmetrizer."""
        from metatomic_ase import SymmetrizedCalculator

        calculator = getattr(model, "calculator", model)
        if not isinstance(calculator, SymmetrizedCalculator):
            calculator = SymmetrizedCalculator(
                calculator, batch_size=settings.batch_size, include_inversion=False
            )
        return PetBatchedModel(calculator)
