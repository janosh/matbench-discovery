import numpy as np
from pymatgen.core import Structure

__author__ = "Janosh Riebesell"
__date__ = "2022-12-02"

np.random.seed(0)  # ensure reproducible structure perturbations


def perturb_structure(struct: Structure, gamma: float = 1.5) -> Structure:
    """Perturb the atomic coordinates of a pymatgen structure

    Args:
        struct (Structure): pymatgen structure to be perturbed

    Returns:
        Structure: Perturbed structure
    """
    perturbed = struct.copy()
    for site in perturbed:
        magnitude = np.random.weibull(gamma)
        vec = np.random.randn(3)  # TODO maybe make func recursive to deal with 0-vector
        vec /= np.linalg.norm(vec)  # unit vector
        site.coords += vec * magnitude
        site.to_unit_cell(in_place=True)

    return perturbed
