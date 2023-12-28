from __future__ import annotations

from typing import TYPE_CHECKING

from matbench_discovery.structure import perturb_structure

if TYPE_CHECKING:
    from pymatgen.core import Structure


def test_perturb_structure(dummy_struct: Structure) -> None:
    perturbed = perturb_structure(dummy_struct)
    assert len(perturbed) == len(dummy_struct)

    for site, new in zip(dummy_struct, perturbed):
        assert site.specie == new.specie
        assert tuple(site.coords) != tuple(new.coords)

    # but different on subsequent calls
    assert perturb_structure(dummy_struct) != perturb_structure(dummy_struct)
