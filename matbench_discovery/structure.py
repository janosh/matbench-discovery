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


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    gamma = 1.5
    samples = np.array([np.random.weibull(gamma) for _ in range(10000)])
    mean = samples.mean()

    # reproduces the dist in https://www.nature.com/articles/s41524-022-00891-8#Fig5
    ax = plt.hist(samples, bins=100)
    # add vertical line at the mean
    plt.axvline(mean, color="gray", linestyle="dashed", linewidth=1)
    # annotate the mean line
    plt.annotate(
        f"{mean = :.2f}",
        xy=(mean, 1),
        # use ax coords for y
        xycoords=("data", "axes fraction"),
        # add text offset
        xytext=(10, -20),
        textcoords="offset points",
    )
