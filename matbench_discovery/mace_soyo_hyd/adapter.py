"""Expose the frozen dimer batching and scalar periodic potential to shared runners."""

import hashlib
from collections.abc import Sequence
from pathlib import Path
from typing import ClassVar

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from .calculator import CanonicalBatchedDimerCalculator, isolated_displacement, key


class HYDCalculator(CanonicalBatchedDimerCalculator):
    """Batch isolated-dimer energy/forces; delegate stress and crystals to scalar ASE.

    The prefetched curve is a fixed input grid, not a reference-data lookup.
    Distances outside that grid are calculated and appended to the same batch.
    """

    implemented_properties: ClassVar[list[str]] = [
        "energy",
        "free_energy",
        "forces",
        "stress",
    ]

    def __init__(self, package_path: str | Path, device: str) -> None:
        """Load the frozen potential with the evaluated fixed-grid batch plan."""
        from matbench_discovery.diatomics import generate_diatomics

        distances = np.geomspace(0.1, 6.0, 119)

        def plan(atoms: Atoms) -> list[Atoms]:
            first, second = atoms.get_chemical_symbols()
            return generate_diatomics(first, second, distances)

        correction = Path(__file__).with_name("frozen-candidate.json")
        super().__init__(package_path, correction, device, plan=plan)

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy", "forces"),
        system_changes: Sequence[str] = all_changes,
    ) -> None:
        """Use canonical dimer batches or the frozen scalar periodic/stress path."""
        atoms = self.atoms if atoms is None else atoms
        if atoms is None:
            raise ValueError("An ASE Atoms object is required")
        target = key(atoms)
        cached = next(iter(self.cached), None)
        off_grid_same_family = (
            cached is not None
            and (target[0], target[2], target[3]) == (cached[0], cached[2], cached[3])
            and target not in self.cached
        )
        if (
            "stress" in properties
            or isolated_displacement(atoms, self.model.cutoff) is None
            or off_grid_same_family
        ):
            Calculator.calculate(self, atoms, properties, system_changes)
            self.model.calculate(self.atoms, properties, system_changes)
            self.results = {
                name: value.copy() if hasattr(value, "copy") else value
                for name, value in self.model.results.items()
            }
            return
        super().calculate(atoms, properties, system_changes)


def make_calculator(device: str, checkpoint: str | Path) -> HYDCalculator:
    """Use the upstream exporter for raw checkpoints, or accept a frozen AOTI file."""
    if not str(checkpoint).endswith(".pt2"):
        from matbench_discovery.calculators import _mace_soyo

        expected = "180784bad2fdee6aec58dae2cdb0d581ad1c4f3d2272d2857bce6315ea43b82b"
        if hashlib.sha256(Path(checkpoint).read_bytes()).hexdigest() != expected:
            raise ValueError("Expected the frozen upstream MACE-SOYO-OAM checkpoint")
        upstream = _mace_soyo("salexmp")(device, checkpoint)
        checkpoint = upstream.package_path
        del upstream
    return HYDCalculator(checkpoint, device)
