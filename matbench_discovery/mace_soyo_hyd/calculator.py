"""Frozen radial correction and canonical-input batched ASE calculator.

Scientific method bodies are copied from the evaluated sources. The source
archive and source-method-verification.json retain their exact provenance.
"""

import json
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import ClassVar

import numpy as np
import torch
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.geometry import find_mic
from ase.stress import full_3x3_to_voigt_6_stress
from mace_soyo import MACESoyoCalculator

BATCH_SIZE = 32
RADII = {
    "H": 0.31,
    "Li": 1.28,
    "C": 0.76,
    "O": 0.66,
    "F": 0.57,
    "Na": 1.66,
    "Al": 1.21,
    "Si": 1.11,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.2,
}
EN = {
    "H": 2.2,
    "Li": 0.98,
    "C": 2.55,
    "O": 3.44,
    "F": 3.98,
    "Na": 0.93,
    "Al": 1.61,
    "Si": 1.9,
    "S": 2.58,
    "Cl": 3.16,
    "Br": 2.96,
}
SCALE_RADII = {
    **RADII,
    "Be": 0.96,
    "B": 0.84,
    "N": 0.71,
    "Mg": 1.41,
    "P": 1.07,
    "K": 2.03,
    "Ca": 1.76,
}
ELECTRONEGATIVITY = {
    **EN,
    "Be": 1.57,
    "B": 2.04,
    "N": 3.04,
    "Mg": 1.31,
    "P": 2.19,
    "K": 0.82,
    "Ca": 1.0,
}
CENTERS = np.array([0.8, 1.0, 1.2, 1.4])
WIDTH = 0.18


def load_coefficients(path: Path) -> np.ndarray:
    """Load and validate the frozen radial parameter array."""
    candidate = json.loads(path.read_text())
    coefficients = np.asarray(candidate["coefficients"], dtype=float)
    if candidate["family"] != "radial" or coefficients.shape != (len(CENTERS) * 4,):
        raise ValueError("Expected the frozen radial candidate")
    if not np.isfinite(coefficients).all():
        raise ValueError("Candidate has non-finite coefficients")
    return coefficients.reshape(len(CENTERS), 4)


def isolated_pair_correction(
    elements: Sequence[str], positions: np.ndarray, coefficients: np.ndarray
) -> tuple[float, np.ndarray]:
    """Return energy and equal-and-opposite forces in eV and eV/Å."""
    a, b = elements
    if any(element not in SCALE_RADII for element in (a, b)):
        return (0.0, np.zeros((2, 3)))
    vector = np.asarray(positions[1]) - np.asarray(positions[0])
    distance = float(np.linalg.norm(vector))
    if not np.isfinite(distance) or distance <= 0.0:
        raise ValueError("Expected two distinct finite atom positions")
    radius = SCALE_RADII[a] + SCALE_RADII[b]
    ena, enb = (ELECTRONEGATIVITY[a], ELECTRONEGATIVITY[b])
    chemistry = np.array(
        [1.0, (ena + enb - 5.0) / 2.0, abs(ena - enb) / 2.0, (radius - 2.0) / 2.0]
    )
    offset = distance / radius - CENTERS
    weighted = np.exp(-0.5 * (offset / WIDTH) ** 2) * (coefficients @ chemistry)
    energy = float(weighted.sum())
    radial = float(weighted @ offset / (radius * WIDTH**2))
    force = radial * vector / distance
    return (energy, np.stack((-force, force)))


def isolated_displacement(atoms: Atoms, cutoff: float) -> np.ndarray | None:
    """Return the pair vector only when the system is an isolated dimer."""
    if len(atoms) != 2:
        return None
    vector = atoms.positions[1] - atoms.positions[0]
    if not atoms.pbc.any():
        return vector
    active_cell = atoms.cell.array[atoms.pbc]
    shortest_bound = np.linalg.svd(active_cell, compute_uv=False)[-1]
    if shortest_bound <= 2 * cutoff:
        return None
    nearest, distance = find_mic(vector, atoms.cell, pbc=atoms.pbc)
    if shortest_bound - distance <= cutoff:
        return None
    return nearest


class MaceSoyoHydCalculator(MACESoyoCalculator):
    """Use the base potential everywhere, with an additive term for isolated dimers."""

    def __init__(
        self,
        package_path: str | Path,
        correction_path: str | Path,
        device: str = "cuda:0",
    ) -> None:
        """Load the frozen potential and radial coefficients."""
        super().__init__(
            package_path=str(package_path), head="salexmp", device=device, use_d3=False
        )
        self.coefficients = load_coefficients(Path(correction_path))

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy", "forces", "stress"),
        system_changes: Sequence[str] = all_changes,
    ) -> None:
        """Populate ASE results for the requested geometry."""
        super().calculate(atoms, properties, system_changes)
        vector = isolated_displacement(self.atoms, float(self.cutoff))
        if vector is None:
            return
        delta_energy, delta_force = isolated_pair_correction(
            self.atoms.get_chemical_symbols(),
            np.stack((np.zeros(3), vector)),
            self.coefficients,
        )
        self.results["energy"] += delta_energy
        if "free_energy" in self.results:
            self.results["free_energy"] += delta_energy
        self.results["forces"] += delta_force
        if self.atoms.pbc.all() and "stress" in self.results:
            stress = -np.outer(vector, delta_force[1]) / self.atoms.get_volume()
            self.results["stress"] += full_3x3_to_voigt_6_stress(stress)


def key(atoms: Atoms) -> tuple[bytes, bytes, bytes, bytes]:
    """Identify an exact ASE input geometry for the prediction cache."""
    return (
        atoms.numbers.tobytes(),
        atoms.positions.tobytes(),
        atoms.cell.array.tobytes(),
        atoms.pbc.tobytes(),
    )


class CanonicalBatchedDimerCalculator(Calculator):
    """Batch compiled calls while retaining the scalar ASE input conversion."""

    implemented_properties: ClassVar[list[str]] = ["energy", "free_energy", "forces"]

    def __init__(
        self,
        package_path: str | Path,
        correction_path: str | Path,
        device: str,
        plan: Callable[[Atoms], list[Atoms]] | None = None,
    ) -> None:
        """Load the frozen potential and radial coefficients."""
        Calculator.__init__(self)
        self.model = MaceSoyoHydCalculator(package_path, correction_path, device)
        self.coefficients = self.model.coefficients
        self.plan = plan
        self.cached = {}

    def predict_base(self, group: list[Atoms]) -> tuple[np.ndarray, np.ndarray]:
        """Concatenate scalar ASE graph inputs without changing their conversion."""
        if self.model._conditioned:  # noqa: SLF001
            raise ValueError(
                "Canonical batching currently requires an unconditioned package"
            )
        inputs = [self.model._atoms_to_inputs(atoms) for atoms in group]  # noqa: SLF001
        z = torch.cat([item[0] for item in inputs])
        fractional = torch.cat([item[1] for item in inputs])
        cells = torch.stack([item[2] for item in inputs])
        batches = torch.cat([item[3] + index for index, item in enumerate(inputs)])
        edges = torch.cat(
            [item[4] + 2 * index for index, item in enumerate(inputs)], dim=1
        )
        images = torch.cat([item[5] for item in inputs])
        args = (z, fractional, cells, batches, edges, images)
        if self.model._dynamic_head:  # noqa: SLF001
            heads = self.model._head_tensor.expand(len(group)).contiguous()  # noqa: SLF001
            pbc = torch.tensor(
                [atoms.pbc.tolist() for atoms in group],
                dtype=torch.bool,
                device=self.model.device,
            )
            args += (heads, pbc)
        with torch.cuda.device(self.model.device), torch.inference_mode():
            energy, force, _ = self.model.aoti_model(*args)
        return (
            energy.cpu().numpy().reshape(-1),
            force.cpu().numpy().reshape(len(group), 2, 3),
        )

    def predict(self, systems: list[Atoms]) -> list[dict[str, float | np.ndarray]]:
        """Evaluate input systems in bounded canonical batches."""
        records = []
        for offset in range(0, len(systems), BATCH_SIZE):
            group = systems[offset : offset + BATCH_SIZE]
            if any(len(atoms) != 2 for atoms in group):
                raise ValueError("Dimer batch requires exactly two atoms per structure")
            energies, forces = self.predict_base(group)
            for atoms, energy, force in zip(group, energies, forces, strict=True):
                vector = isolated_displacement(atoms, self.model.cutoff)
                if vector is not None:
                    delta_energy, delta_force = isolated_pair_correction(
                        atoms.get_chemical_symbols(),
                        np.stack((np.zeros(3), vector)),
                        self.coefficients,
                    )
                    energy += delta_energy  # noqa: PLW2901 - preserve frozen arithmetic
                    force = force.copy()  # noqa: PLW2901 - preserve frozen arithmetic
                    force += delta_force  # noqa: PLW2901 - preserve frozen arithmetic
                records.append(
                    {
                        "energy": float(energy),
                        "free_energy": float(energy),
                        "forces": force.copy(),
                    }
                )
        return records

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy", "forces"),
        system_changes: Sequence[str] = all_changes,
    ) -> None:
        """Populate ASE results for the requested geometry."""
        super().calculate(atoms, properties, system_changes)
        target = key(self.atoms)
        if target not in self.cached:
            group = self.plan(self.atoms) if self.plan else [self.atoms]
            keys = [key(atoms) for atoms in group]
            if target not in keys:
                group.append(self.atoms)
                keys.append(target)
            self.cached = dict(zip(keys, self.predict(group), strict=True))
        self.results = self.cached[target]
