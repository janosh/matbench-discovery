from __future__ import annotations

import os
import urllib.error
from collections.abc import Sequence
from glob import glob
from pathlib import Path
from typing import Any, Callable

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from tqdm import tqdm

from matbench_discovery import ROOT

# repo URL to raw files on GitHub
RAW_REPO_URL = "https://raw.githubusercontent.com/janosh/matbench-discovery"
# directory to cache downloaded data files
default_cache_dir = os.path.expanduser("~/.cache/matbench-discovery")


class Files(dict):  # type: ignore
    """Files instance inherits from dict so that .values(), items(), etc. are supported
    but also allows accessing attributes by dot notation. E.g. FILES.wbm_summary instead
    of FILES["wbm_summary"]. This enables tab completion in IDEs and auto-updating
    attribute names across the code base when changing the name of an attribute. Every
    subclass must set the _root attribute to a path that serves as the root directory
    w.r.t. which all files will be turned into absolute paths. The _key_map attribute
    can be used to map attribute names to different names in the dict. Useful if you
    want to have keys like 'foo+bar' that are not valid Python identifiers.
    """

    def __init__(self) -> None:
        """Create a Files instance."""
        key_map = getattr(self, "_key_map", {})
        dct = {
            key_map.get(key, key): f"{self._root}{file}"  # type: ignore
            for key, file in type(self).__dict__.items()
            if not key.startswith("_")
        }
        self.__dict__ = dct
        super().__init__(dct)


class DataFiles(Files):
    """Data files provided by Matbench Discovery.
    See https://janosh.github.io/matbench-discovery/contribute for data descriptions.
    """

    _root = f"{ROOT}/data/"

    mp_computed_structure_entries = (
        "mp/2023-02-07-mp-computed-structure-entries.json.gz"
    )
    mp_elemental_ref_entries = "mp/2022-09-19-mp-elemental-reference-entries.json"
    mp_energies = "mp/2023-01-10-mp-energies.csv"
    mp_patched_phase_diagram = "mp/2023-02-07-ppd-mp.pkl.gz"
    wbm_computed_structure_entries = (
        "wbm/2022-10-19-wbm-computed-structure-entries.json.bz2"
    )
    wbm_initial_structures = "wbm/2022-10-19-wbm-init-structs.json.bz2"
    wbm_cses_plus_init_structs = (
        "wbm/2022-10-19-wbm-computed-structure-entries+init-structs.json.bz2"
    )
    wbm_summary = "wbm/2022-10-19-wbm-summary.csv"


DATA_FILES = DataFiles()


def as_dict_handler(obj: Any) -> dict[str, Any] | None:
    """Pass this to json.dump(default=) or as pandas.to_json(default_handler=) to
    serialize Python classes with as_dict(). Warning: Objects without a as_dict() method
    are replaced with None in the serialized data.
    """
    try:
        return obj.as_dict()  # all MSONable objects implement as_dict()
    except AttributeError:
        return None  # replace unhandled objects with None in serialized data
        # removes e.g. non-serializable AseAtoms from M3GNet relaxation trajectories


def load_train_test(
    data_names: str | Sequence[str] = ("summary",),
    version: str = "1.0.0",
    cache_dir: str | Path = default_cache_dir,
    hydrate: bool = False,
    **kwargs: Any,
) -> pd.DataFrame:
    """Download parts of or the full MP training data and WBM test data as pandas
    DataFrames. The full training and test sets are each about ~500 MB as compressed
    JSON which will be cached locally to cache_dir for faster re-loading unless
    cache_dir is set to None.

    See matbench_discovery.data.DATA_FILES for recognized data keys. See
    https://janosh.github.io/matbench-discovery/contribute for descriptions.

    Args:
        data_names (str | list[str], optional): Which parts of the MP/WBM data to load.
            Can be any subset of the above data names or 'all'. Defaults to ["summary"].
        version (str, optional): Which version of the dataset to load. Defaults to
            '1.0.0'. Can be any git tag, branch or commit hash.
        cache_dir (str, optional): Where to cache data files on local drive. Defaults to
            '~/.cache/matbench-discovery'. Set to None to disable caching.
        hydrate (bool, optional): Whether to hydrate pymatgen objects. If False,
            Structures and ComputedStructureEntries are returned as dictionaries which
            can be hydrated on-demand with df.col.map(Structure.from_dict). Defaults to
            False as it noticeably increases load time.
        **kwargs: Additional keyword arguments passed to pandas.read_json or read_csv,
            depending on which file is loaded.

    Raises:
        ValueError: On bad version number or bad data names.

    Returns:
        pd.DataFrame: Single dataframe or dictionary of dfs if
        multiple data were requested.
    """
    if data_names == "all":
        data_names = list(DATA_FILES)
    elif isinstance(data_names, str):
        data_names = [data_names]

    if missing := set(data_names) - set(DATA_FILES):
        raise ValueError(f"{missing} must be subset of {set(DATA_FILES)}")

    dfs = {}
    for key in data_names:
        file = DataFiles.__dict__[key]
        reader = pd.read_csv if file.endswith(".csv") else pd.read_json

        cache_path = f"{cache_dir}/{version}/{file}"
        if os.path.isfile(cache_path):
            print(f"Loading {key!r} from cached file at {cache_path!r}")
            df = reader(cache_path, **kwargs)
        else:
            url = f"{RAW_REPO_URL}/{version}/data/{file}"
            print(f"Downloading {key!r} from {url}")
            try:
                df = reader(url)
            except urllib.error.HTTPError as exc:
                raise ValueError(f"Bad {url=}") from exc
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

    if len(data_names) == 1:
        return dfs[data_names[0]]
    return dfs


def glob_to_df(
    pattern: str,
    reader: Callable[[Any], pd.DataFrame] = None,
    pbar: bool = True,
    **kwargs: Any,
) -> pd.DataFrame:
    """Combine data files matching a glob pattern into a single dataframe.

    Args:
        pattern (str): Glob file pattern.
        reader (Callable[[Any], pd.DataFrame], optional): Function that loads data from
            disk. Defaults to pd.read_csv if ".csv" in pattern else pd.read_json.
        pbar (bool, optional): Whether to show progress bar. Defaults to True.
        **kwargs: Keyword arguments passed to reader (i.e. pd.read_csv or pd.read_json).

    Returns:
        pd.DataFrame: Combined dataframe.
    """
    reader = reader or pd.read_csv if ".csv" in pattern else pd.read_json

    # prefix pattern with ROOT if not absolute path
    files = glob(pattern)
    if len(files) == 0:
        raise FileNotFoundError(f"No files matching glob {pattern=}")

    sub_dfs = {}  # used to join slurm job array results into single df
    for file in tqdm(files, disable=not pbar):
        df = reader(file, **kwargs)
        sub_dfs[file] = df

    return pd.concat(sub_dfs.values())


df_wbm = load_train_test("wbm_summary")
df_wbm.index = df_wbm.material_id
