from pymatgen.core import Structure

from matbench_discovery.structure import perturb_structure


def test_perturb_structure(dummy_struct: Structure) -> None:
    perturbed = perturb_structure(dummy_struct)
    assert len(perturbed) == len(dummy_struct)

    for site, new in zip(dummy_struct, perturbed, strict=True):
        assert site.specie == new.specie
        assert tuple(site.coords) != tuple(new.coords)

    # but different on subsequent calls
    assert perturb_structure(dummy_struct) != perturb_structure(dummy_struct)
