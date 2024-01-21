from __future__ import annotations

import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry

from matbench_discovery import Key


@pytest.fixture()
def dummy_df_serialized(dummy_struct: Structure) -> pd.DataFrame:
    # create a dummy df with a structure column on which to test (de-)serialization
    df = pd.DataFrame(dict(material_id=range(5), structure=[dummy_struct] * 5))
    df[Key.volume] = [x.volume for x in df.structure]
    df[Key.struct] = [x.as_dict() for x in df.structure]
    cse_dict = ComputedStructureEntry(dummy_struct, 0).as_dict()
    df[Key.cse] = [cse_dict] * len(df)
    return df


@pytest.fixture()
def dummy_struct() -> Structure:
    return Structure(
        lattice=Lattice.cubic(4.2),
        species=("Fe", "O"),
        coords=((0, 0, 0), (0.5, 0.5, 0.5)),
    )
