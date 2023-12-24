from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from matbench_discovery.structure import perturb_structure

if TYPE_CHECKING:
    from pymatgen.core import Structure


def test_perturb_structure(dummy_struct: Structure) -> None:
    np.random.seed(0)
    perturbed = perturb_structure(dummy_struct)
    assert len(perturbed) == len(dummy_struct)

    for site, new in zip(dummy_struct, perturbed):
        assert site.specie == new.specie
        assert tuple(site.coords) != tuple(new.coords)

    # test that the perturbation is reproducible
    np.random.seed(0)
    assert perturbed == perturb_structure(dummy_struct)
    # but different on subsequent calls
    assert perturb_structure(dummy_struct) != perturb_structure(dummy_struct)
