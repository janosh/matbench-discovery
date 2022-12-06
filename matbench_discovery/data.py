from __future__ import annotations

import os
from collections.abc import Generator, Sequence
from typing import Any

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from tqdm import tqdm

DATA_FILENAMES = {
    "wbm-summary": "wbm/2022-10-19-wbm-summary.csv",
    "wbm-initial-structures": "wbm/2022-10-19-wbm-init-structs.json.bz2",
    "wbm-computed-structure-entries": "wbm/2022-10-19-wbm-cses.json.bz2",
    "mp-energies": "mp/2022-08-13-mp-energies.json.gz",
    "mp-computed-structure-entries": "mp/2022-09-16-mp-computed-structure-entries.json.gz",
    "mp-patched-phase-diagram": "mp/2022-09-18-ppd-mp.pkl.gz",
    "mp-elemental-ref-energies": "mp/2022-09-19-mp-elemental-ref-energies.json",
}

RAW_REPO_URL = "https://raw.githubusercontent.com/janosh/matbench-discovery"
default_cache_dir = os.path.expanduser("~/.cache/matbench-discovery")


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


def load_train_test(
    parts: str | Sequence[str] = ("summary",),
    version: int = 1,
    cache_dir: str | None = default_cache_dir,
    hydrate: bool = False,
) -> pd.DataFrame | dict[str, pd.DataFrame]:
    """Download the MP training data and WBM test data in parts or in full as pandas
    DataFrames. The full training and test sets are each about ~500 MB as compressed
    JSON will be cached locally for faster re-loading unless cache_dir is set to None.

    Hint: Import DATA_FILES from the same module as this function and
    print(list(DATA_FILES)) to see permissible data names.

    Args:
        parts (str | list[str], optional): Which parts of the MP/WBM dataset to load.
            Can be any subset of list(DATA_FILES). Defaults to ["summary"], a dataframe
            with columns for material properties like VASP energy, formation energy,
            energy above the convex hull (3 columns with old, new and no Materials
            Project energy corrections applied for each), volume, band gap, number of
            sites per unit cell, and more.
        version (int, optional): Which version of the dataset to load. Defaults to 1
            (currently the only available option).
        cache_dir (str, optional): Where to cache data files on local drive. Defaults to
            '~/.cache/matbench-discovery'. Set to None to disable caching.
        hydrate (bool, optional): Whether to hydrate pymatgen objects. If False,
            Structures and ComputedStructureEntries are returned as dictionaries which
            can be hydrated on-demand with df.col.map(Structure.from_dict). Defaults to
            False as it noticeably increases load time.

    Raises:
        ValueError: On bad version number or bad part names.

    Returns:
        pd.DataFrame | dict[str, pd.DataFrame]: Single dataframe of dictionary of
        multiple data parts were requested.
    """
    if parts == "all":
        parts = list(DATA_FILENAMES)
    elif isinstance(parts, str):
        parts = [parts]

    if version != 1:
        raise ValueError(f"Only version 1 currently available, got {version=}")
    if missing := set(parts) - set(DATA_FILENAMES):
        raise ValueError(f"{missing} must be subset of {set(DATA_FILENAMES)}")

    dfs = {}
    for key in parts:
        file = DATA_FILENAMES[key]
        reader = pd.read_csv if file.endswith(".csv") else pd.read_json

        cache_path = f"{cache_dir}/{file}"
        if os.path.isfile(cache_path):
            df = reader(cache_path)
        else:
            url = f"{RAW_REPO_URL}/{version}.0.0/data/{file}"
            print(f"Downloading {key} from {url}")
            df = reader(url)
            if cache_dir and not os.path.isfile(cache_path):
                os.makedirs(os.path.dirname(cache_path), exist_ok=True)
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
