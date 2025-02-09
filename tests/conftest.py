import numpy as np
import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key


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
