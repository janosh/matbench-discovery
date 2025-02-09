import pytest
from pymatgen.core import Structure

from matbench_discovery.structure import perturb_structure


@pytest.mark.parametrize("gamma", [0.5, 1.0, 1.5, 2.0])
def test_perturb_structure(dummy_struct: Structure, gamma: float) -> None:
    perturbed = perturb_structure(dummy_struct, gamma=gamma)
    assert len(perturbed) == len(dummy_struct)

    for site, new in zip(dummy_struct, perturbed, strict=True):
        assert site.specie == new.specie
        assert tuple(site.coords) != tuple(new.coords)

    # Check that different calls produce different results
    assert perturb_structure(dummy_struct, gamma=gamma) != perturb_structure(
        dummy_struct, gamma=gamma
    )
