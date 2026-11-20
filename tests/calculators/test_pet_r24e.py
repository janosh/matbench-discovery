"""Coordinate and energy-reference contracts for the R24E calculator."""

from collections.abc import Sequence

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.stress import full_3x3_to_voigt_6_stress

from matbench_discovery.pet_r24e import R24E, rotations


class Anisotropic(Calculator):
    """A direction-dependent harmonic energy with analytic force and stress."""

    implemented_properties = ("energy", "forces", "stress")

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy",),
        system_changes: Sequence[str] = all_changes,
    ) -> None:
        """Compute the requested ASE properties in the input frame."""
        super().calculate(atoms, properties, system_changes)
        pos = self.atoms.positions
        forces = -pos * np.array([1.0, 2.0, 6.0])
        self.results = {"energy": -0.5 * (pos * forces).sum(), "forces": forces}
        if "stress" in properties:
            self.results["stress"] = full_3x3_to_voigt_6_stress(
                -pos.T @ forces / self.atoms.get_volume()
            )


@pytest.mark.parametrize("periodic", [False, True])
def test_orientation_and_element_reference(periodic: bool) -> None:
    atoms = Atoms(
        "HSi",
        positions=[[0.2, 0.1, 0.4], [1.3, 0.5, 0.8]],
        cell=np.eye(3) * 5 if periodic else None,
        pbc=periodic,
    )
    z = np.zeros(119)
    z[[1, 14]] = [-1.25, 0.75]
    atoms.calc = R24E(Anisotropic(), intercept=0.1, offsets=z)
    assert atoms.get_potential_energy() == pytest.approx(
        1.5 * (atoms.positions**2).sum() - 0.3
    )
    np.testing.assert_allclose(atoms.get_forces(), -3 * atoms.positions, atol=1e-11)
    if periodic:
        expected = full_3x3_to_voigt_6_stress(
            3 * atoms.positions.T @ atoms.positions / atoms.get_volume()
        )
        np.testing.assert_allclose(atoms.get_stress(), expected, atol=1e-11)


def test_rotations_and_registry() -> None:
    from matbench_discovery.calculators import CALCULATORS

    matrices = rotations()
    assert set(np.round(np.linalg.det(matrices), 8)) == {-1, 1}
    assert len({tuple(r.ravel()) for r in matrices}) == 24
    assert CALCULATORS["pet_omat_xs_r24e"].auto_checkpoint
