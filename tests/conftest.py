from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry

from matbench_discovery import Key

rng = np.random.default_rng(0)


@pytest.fixture()
def dummy_struct() -> Structure:
    return Structure(
        lattice=Lattice.cubic(4.2),
        species=("Fe", "O"),
        coords=((0, 0, 0), (0.5, 0.5, 0.5)),
    )


@pytest.fixture()
def df_float() -> pd.DataFrame:
    return pd.DataFrame(rng.random(size=(10, 5)), columns=[*"ABCDE"])


@pytest.fixture()
def df_mixed() -> pd.DataFrame:
    return pd.DataFrame(
        dict(
            A=rng.random(size=10),
            B=rng.choice([True, False], size=10),
            C=rng.choice([*"abcdef"], size=10),
        )
    )


@pytest.fixture()
def df_with_pmg_objects(dummy_struct: Structure) -> pd.DataFrame:
    # create a dummy df with a structure column on which to test (de-)serialization
    df = pd.DataFrame(dict(material_id=range(5), structure=[dummy_struct] * 5))
    df[Key.volume] = [x.volume for x in df.structure]
    df[Key.struct] = [x.as_dict() for x in df.structure]
    cse_dict = ComputedStructureEntry(dummy_struct, 0).as_dict()
    df[Key.cse] = [cse_dict] * len(df)
    return df
