import numpy as np
import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key

from matbench_discovery.metrics.diatomics import DiatomicCurves


@pytest.fixture
def dummy_struct() -> Structure:
    return Structure(
        lattice=Lattice.cubic(4.2),
        species=("Fe", "O"),
        coords=((0, 0, 0), (0.5, 0.5, 0.5)),
    )


@pytest.fixture
def cubic_struct() -> Structure:
    lattice = Lattice.cubic(4.2)
    return Structure(lattice, ["Si", "Si"], [[0, 0, 0], [0.5, 0.5, 0.5]])


@pytest.fixture
def tetragonal_struct() -> Structure:
    lattice = Lattice.tetragonal(a=4, c=6)
    coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5]]
    return Structure(lattice, ["Ti", "O", "O"], coords)


@pytest.fixture
def monoclinic_struct() -> Structure:
    lattice = Lattice.monoclinic(a=5, b=4, c=6, beta=100)
    x1 = [0, 0, 0]
    x2 = [0.25, 0.25, 0.25]
    x3 = [0.75, 0.75, 0.25]
    x4 = [0.25, 0.75, 0.75]
    x5 = [0.75, 0.25, 0.75]
    return Structure(lattice, ["P", "O", "O", "O", "O"], [x1, x2, x3, x4, x5])


@pytest.fixture
def df_float() -> pd.DataFrame:
    rng = np.random.default_rng(seed=0)

    return pd.DataFrame(rng.normal(size=(10, 5)), columns=[*"ABCDE"])


@pytest.fixture
def df_mixed() -> pd.DataFrame:
    rng = np.random.default_rng(seed=0)

    floats = rng.random(size=10)
    bools = rng.choice([True, False], size=10)
    strings = rng.choice([*"abcdef"], size=10)
    return pd.DataFrame(dict(floats=floats, bools=bools, strings=strings))


@pytest.fixture
def df_with_pmg_objects(dummy_struct: Structure) -> pd.DataFrame:
    # create a dummy df with a structure column on which to test (de-)serialization
    df_dummy = pd.DataFrame(dict(material_id=range(5), structure=[dummy_struct] * 5))
    df_dummy[Key.volume] = [x.volume for x in df_dummy.structure]
    df_dummy[Key.structure] = [x.as_dict() for x in df_dummy.structure]
    cse_dict = ComputedStructureEntry(dummy_struct, 0).as_dict()
    df_dummy[Key.computed_structure_entry] = [cse_dict] * len(df_dummy)
    return df_dummy


@pytest.fixture
def pred_ref_diatomic_curves() -> tuple[DiatomicCurves, DiatomicCurves]:
    """Create test data for diatomic curves.

    Returns:
        TestCurves: Reference and predicted curves for each element.
    """
    # Simple test case: H2 molecule with slightly different curves
    distances = np.linspace(0.5, 5.0, 100)
    # Reference: Simple Morse potential
    e_h_ref = 5 * (1 - np.exp(-2 * (distances - 1.5))) ** 2 - 5
    e_he_ref = 4 * (1 - np.exp(-2 * (distances - 1.5))) ** 2 - 4
    # Prediction: Slightly perturbed Morse potential
    e_h_pred = 5.2 * (1 - np.exp(-2.1 * (distances - 1.48))) ** 2 - 5.1
    e_he_pred = 4.2 * (1 - np.exp(-2.1 * (distances - 1.48))) ** 2 - 4.1

    np_rng = np.random.default_rng(seed=0)

    # forces have shape (2, 3), 2 atoms with 3 force components
    f_h_ref = 0.5 * np_rng.normal(size=(len(distances), 2, 3))
    f_h_pred = 0.55 * np_rng.normal(size=(len(distances), 2, 3))
    f_he_ref = 0.1 * np_rng.normal(size=(len(distances), 2, 3))
    f_he_pred = 0.15 * np_rng.normal(size=(len(distances), 2, 3))

    ref_curves = {
        "distances": distances,
        "homo-nuclear": {
            "H": {"energies": e_h_ref, "forces": f_h_ref},
            "He": {"energies": e_he_ref, "forces": f_he_ref},
        },
    }
    pred_curves = {
        "distances": distances,
        "homo-nuclear": {
            "H": {"energies": e_h_pred, "forces": f_h_pred},
            "He": {"energies": e_he_pred, "forces": f_he_pred},
        },
    }
    return DiatomicCurves.from_dict(ref_curves), DiatomicCurves.from_dict(pred_curves)
