from __future__ import annotations

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from matbench_discovery.structure import perturb_structure


@pytest.fixture
def struct() -> Structure:
    return Structure(
        lattice=Lattice.cubic(5),
        species=("Fe", "O"),
        coords=((0, 0, 0), (0.5, 0.5, 0.5)),
    )


def test_perturb_structure(struct: Structure) -> None:
    np.random.seed(0)
    perturbed = perturb_structure(struct)
    assert len(perturbed) == len(struct)

    for site, new in zip(struct, perturbed):
        assert site.specie == new.specie
        assert tuple(site.coords) != tuple(new.coords)

    # test that the perturbation is reproducible
    np.random.seed(0)
    assert perturbed == perturb_structure(struct)
    # but different on subsequent calls
    assert perturb_structure(struct) != perturb_structure(struct)
