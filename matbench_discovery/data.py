from __future__ import annotations

import os
from collections.abc import Generator, Sequence
from typing import Any

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from tqdm import tqdm

data_files = {
    "summary": "2022-10-19-wbm-summary.csv",
    "initial-structures": "2022-10-19-wbm-init-structs.json.bz2",
    "computed-structure-entries": "2022-10-19-wbm-cses.json.bz2",
}

base_url = "https://raw.githubusercontent.com/janosh/matbench-discovery/main/data/wbm"
default_cache_loc = os.path.expanduser("~/.cache/matbench-discovery")


def chunks(xs: Sequence[Any], n: int) -> Generator[Sequence[Any], None, None]:
    return (xs[i : i + n] for i in range(0, len(xs), n))


def as_dict_handler(obj: Any) -> dict[str, Any] | None:
    """Pass this to json.dump(default=) or as pandas.to_json(default_handler=) to
    convert Python classes with a as_dict() method to dictionaries on serialization.
    Objects without a as_dict() method are replaced with None in the serialized data.
    """
    try:
        return obj.as_dict()  # all MSONable objects implement as_dict()
    except AttributeError:
        return None  # replace unhandled objects with None in serialized data
        # removes e.g. non-serializable AseAtoms from M3GNet relaxation trajectories


def load_wbm(
    parts: Sequence[str] = ("summary",),
    version: int = 1,
    cache_dir: str | None = default_cache_loc,
    hydrate: bool = False,
) -> pd.DataFrame | dict[str, pd.DataFrame]:
    """_summary_

    Args:
        parts (str, optional): Which parts of the WBM dataset to load. Can be any subset
            of {'summary', 'initial-structures', 'computed-structure-entries'}. Defaults
            to ["summary"], a dataframe with columns for material properties like VASP
            energy, formation energy, energy above the convex hull (3 columns with old,
            new and no Materials Project energy corrections applied for each), volume,
            band gap, number of sites per unit cell, and more.
        version (int, optional): Which version of the dataset to load. Defaults to 1
            (currently the only available option).
        cache_dir (str, optional): Where to cache data files on local drive. Defaults to
            '~/.cache/matbench-discovery'. Set to None to disable caching.
        hydrate (bool, optional): Whether to hydrate pymatgen objects. If False,
            Structures and ComputedStructureEntries are returned as dictionaries which
            can be hydrated on-demand with df.col.map(Structure.from_dict). Defaults to
            False as it noticeably increases load time.

    Raises:
        ValueError: On bad version or bad keys for which data parts to load.

    Returns:
        pd.DataFrame | dict[str, pd.DataFrame]: Single dataframe of dictionary of
        multiple data parts were requested.
    """
    if version != 1:
        raise ValueError(f"Only version 1 currently available, got {version=}")
    if missing := set(parts) - set(data_files):
        raise ValueError(f"{missing} must be subset of {set(data_files)}")

    dfs = {}
    for key in parts:
        file = data_files[key]
        reader = pd.read_csv if file.endswith(".csv") else pd.read_json

        cache_path = f"{cache_dir}/{file}"
        if os.path.isfile(cache_path):
            df = reader(cache_path)
        else:
            url = f"{base_url}/{file}"
            print(f"Downloading {url=}")
            df = reader(url)
            if cache_dir and not os.path.isfile(cache_path):
                os.makedirs(cache_dir, exist_ok=True)
                if ".csv" in file:
                    df.to_csv(cache_path)
                elif ".json" in file:
                    df.reset_index().to_json(
                        cache_path, default_handler=as_dict_handler
                    )
                else:
                    raise ValueError(f"Unexpected file type {file}")

        df = df.set_index("material_id")
        if hydrate:
            for col in df:
                if not isinstance(df[col].iloc[0], dict):
                    continue
                try:
                    df[col] = [
                        ComputedStructureEntry.from_dict(d)
                        for d in tqdm(df[col], desc=col)
                    ]
                except Exception:
                    df[col] = [Structure.from_dict(d) for d in tqdm(df[col], desc=col)]

        dfs[key] = df

    if len(parts) == 1:
        return dfs[parts[0]]
    return dfs
