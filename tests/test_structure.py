from __future__ import annotations

import numpy as np
from pymatgen.core import Structure

from matbench_discovery.structure import perturb_structure


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
