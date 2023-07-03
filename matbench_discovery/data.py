from __future__ import annotations

import gzip
import json
import os
import pickle
import sys
import urllib.error
import urllib.request
from glob import glob
from pathlib import Path
from typing import Any, Callable

import pandas as pd
from monty.json import MontyDecoder
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram
from tqdm import tqdm

from matbench_discovery import FIGSHARE

# repo URL to raw files on GitHub
RAW_REPO_URL = "https://github.com/janosh/matbench-discovery/raw"
figshare_versions = sorted(
    x.split(os.path.sep)[-1].split(".json")[0] for x in glob(f"{FIGSHARE}/*.json")
)
# directory to cache downloaded data files
default_cache_dir = os.getenv(
    "MATBENCH_DISCOVERY_CACHE_DIR",
    f"{os.path.expanduser('~/.cache/matbench-discovery')}/{figshare_versions[-1]}",
)


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


def load(
    data_key: str,
    version: str = figshare_versions[-1],
    cache_dir: str | Path = default_cache_dir,
    hydrate: bool = False,
    **kwargs: Any,
) -> pd.DataFrame | PatchedPhaseDiagram:
    """Download parts of or the full MP training data and WBM test data as pandas
    DataFrames. The full training and test sets are each about ~500 MB as compressed
    JSON which will be cached locally to cache_dir for faster re-loading unless
    cache_dir is set to None.

    See matbench_discovery.data.DATA_FILES for recognized data keys. See [here]
    (https://janosh.github.io/matbench-discovery/contribute#--direct-download) for file
    descriptions.

    Args:
        data_key (str): Which parts of the MP/WBM data to load. Must be one of
            list(DATA_FILES).
        version (str, optional): Which version of the dataset to load. Defaults to
            latest version of data files published to Figshare. Pass any invalid version
            to see valid options.
        cache_dir (str, optional): Where to cache data files on local drive. Defaults to
            '~/.cache/matbench-discovery'. Set to None to disable caching.
        hydrate (bool, optional): Whether to hydrate pymatgen objects. If False,
            Structures and ComputedStructureEntries are returned as dictionaries which
            can be hydrated on-demand with df.col.map(Structure.from_dict). Defaults to
            False as it noticeably increases load time.
        **kwargs: Additional keyword arguments passed to pandas.read_json or read_csv,
            depending on which file is loaded.

    Raises:
        ValueError: On bad version number or bad data_key.

    Returns:
        pd.DataFrame: Single dataframe or dictionary of dfs if multiple data requested.
    """
    if version not in figshare_versions:
        raise ValueError(f"Unexpected {version=}. Must be one of {figshare_versions}.")

    if not isinstance(data_key, str) or data_key not in DATA_FILES:
        raise ValueError(f"Unknown {data_key=}, must be one of {list(DATA_FILES)}.")

    with open(f"{FIGSHARE}/{version}.json") as json_file:
        file_urls = json.load(json_file)

    file = DataFiles.__dict__[data_key]

    cache_path = f"{cache_dir}/{file}"
    if not os.path.isfile(cache_path):  # download from Figshare URL
        url = file_urls[data_key][0]
        print(f"Downloading {data_key!r} from {url}")
        try:
            # ensure directory exists
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            # download and save to disk
            urllib.request.urlretrieve(url, cache_path)
            print(f"Cached {data_key!r} to {cache_path!r}")
        except urllib.error.HTTPError as exc:
            raise ValueError(f"Bad {url=}") from exc
        except Exception:
            print(f"\n\nvariable dump:\n{file=},\n{url=}")
            raise

    print(f"Loading {data_key!r} from cached file at {cache_path!r}")
    if ".pkl" in file:  # handle key='mp_patched_phase_diagram' separately
        with gzip.open(cache_path, "rb") as zip_file:
            return pickle.load(zip_file)

    csv_ext = (".csv", ".csv.gz", ".csv.bz2")
    reader = pd.read_csv if file.endswith(csv_ext) else pd.read_json
    try:
        df = reader(cache_path, **kwargs)
    except Exception:
        print(f"\n\nvariable dump:\n{file=},\n{reader=}\n{kwargs=}")
        raise

    if "material_id" in df:
        df = df.set_index("material_id")
    if hydrate:
        for col in df:
            if not isinstance(df[col].iloc[0], dict):
                continue
            try:
                # convert dicts to pymatgen Structures and ComputedStructureEntrys
                df[col] = [
                    MontyDecoder().process_decoded(dct)
                    for dct in tqdm(df[col], desc=col)
                ]
            except Exception:
                print(f"\n\nvariable dump:\n{col=},\n{df[col]=}")
                raise

    return df


def glob_to_df(
    pattern: str,
    reader: Callable[[Any], pd.DataFrame] | None = None,
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
    reader = reader or pd.read_csv if ".csv" in pattern.lower() else pd.read_json

    # prefix pattern with ROOT if not absolute path
    files = glob(pattern)
    if len(files) == 0:
        raise FileNotFoundError(f"No files matching glob {pattern=}")

    sub_dfs = {}  # used to join slurm job array results into single df
    for file in tqdm(files, disable=not pbar):
        df = reader(file, **kwargs)
        sub_dfs[file] = df

    return pd.concat(sub_dfs.values())


class Files(dict):  # type: ignore
    """Files instance inherits from dict so that .values(), items(), etc. are supported
    but also allows accessing attributes by dot notation. E.g. FILES.wbm_summary instead
    of FILES["wbm_summary"]. This enables tab completion in IDEs and auto-updating
    attribute names across the code base when changing the key of a file.
    """

    def __init__(
        self, root: str = default_cache_dir, key_map: dict[str, str] | None = None
    ) -> None:
        """Create a Files instance.

        Args:
            root (str, optional): Root directory used to absolufy every file path.
            Defaults to '~/.cache/matbench-discovery/[latest_figshare_release]' where
                latest_figshare_release is e.g. 1.0.0. Can also be set through env var
                MATBENCH_DISCOVERY_CACHE_DIR.
            key_map (dict[str, str], optional): Mapping from attribute names to keys in
                the dict. Useful if you want to have keys like 'foo+bar' that are not
                valid Python identifiers. Defaults to None.
        """
        self._on_not_found: Callable[[str, str], None] | None = None
        rel_paths = {
            (key_map or {}).get(key, key): file
            for key, file in type(self).__dict__.items()
            if not key.startswith("_")
        }
        abs_paths = {key: f"{root}/{file}" for key, file in rel_paths.items()}
        self.__dict__ = abs_paths
        super().__init__(abs_paths)

    def __getattribute__(self, key: str) -> str:
        """Override __getattr__ to check if file corresponding to key exists."""
        val = super().__getattribute__(key)
        if key in self and not os.path.isfile(val):
            msg = f"Warning: {val!r} associated with {key=} does not exist."
            if self._on_not_found:
                self._on_not_found(key, msg)
            else:
                print(msg, file=sys.stderr)
        return val


class DataFiles(Files):
    """Data files provided by Matbench Discovery.
    See https://janosh.github.io/matbench-discovery/contribute for data descriptions.
    """

    def _on_not_found(self, key: str, msg: str) -> None:  # type: ignore[override]
        msg += (
            " Would you like to download it now using matbench_discovery."
            f"data.load({key!r}). This will cache the file for future use."
        )

        # default to 'y' if not in interactive session, and user can't answer
        answer = "" if sys.stdin.isatty() else "y"
        while answer not in ("y", "n"):
            answer = input(f"{msg} [y/n] ").lower().strip()
        if answer == "y":
            load(key)  # download and cache data file

    mp_computed_structure_entries = (
        "mp/2023-02-07-mp-computed-structure-entries.json.gz"
    )
    mp_elemental_ref_entries = "mp/2023-02-07-mp-elemental-reference-entries.json.gz"
    mp_energies = "mp/2023-01-10-mp-energies.csv.gz"
    mp_patched_phase_diagram = "mp/2023-02-07-ppd-mp.pkl.gz"
    wbm_computed_structure_entries = (
        "wbm/2022-10-19-wbm-computed-structure-entries.json.bz2"
    )
    wbm_initial_structures = "wbm/2022-10-19-wbm-init-structs.json.bz2"
    wbm_cses_plus_init_structs = (
        "wbm/2022-10-19-wbm-computed-structure-entries+init-structs.json.bz2"
    )
    wbm_summary = "wbm/2022-10-19-wbm-summary.csv.gz"


# data files can be downloaded and cached with matbench_discovery.data.load()
DATA_FILES = DataFiles()


df_wbm = load("wbm_summary")
df_wbm["material_id"] = df_wbm.index
