"""Download, cache and hydrate data files from the Matbench Discovery Figshare article.

https://figshare.com/articles/dataset/22715158
"""

import os
import sys
from collections.abc import Callable
from enum import StrEnum
from glob import glob
from typing import Any, Self

import pandas as pd
import requests
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery import DATA_DIR, FIGSHARE_DIR, pkg_is_editable

# ruff: noqa: T201

# repo URL to raw files on GitHub
RAW_REPO_URL = "https://github.com/janosh/matbench-discovery/raw"
figshare_versions = sorted(
    x.split(os.path.sep)[-1].split(".json")[0] for x in glob(f"{FIGSHARE_DIR}/*.json")
)
# directory to cache downloaded data files
DEFAULT_CACHE_DIR = os.getenv(
    "MATBENCH_DISCOVERY_CACHE_DIR",
    DATA_DIR if pkg_is_editable else os.path.expanduser("~/.cache/matbench-discovery"),
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


def glob_to_df(
    pattern: str,
    *,
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
        df_i = reader(file, **kwargs)
        sub_dfs[file] = df_i

    return pd.concat(sub_dfs.values())


def download_file(file_path: str, url: str) -> None:
    """Download the file from the given URL to the given file path."""
    file_dir = os.path.dirname(file_path)
    os.makedirs(file_dir, exist_ok=True)
    try:
        response = requests.get(url, timeout=5)

        response.raise_for_status()

        with open(file_path, "wb") as file:
            file.write(response.content)
    except requests.exceptions.RequestException as exc:
        print(f"Error downloading {url=}\nto {file_path=}.\n{exc!s}")


class Files(StrEnum):
    """Enum of data files with associated file directories and URLs."""

    def __new__(cls, file_path: str, url: str | None = None) -> Self:
        """Create a new member of the FileUrls enum with a given URL where to load the
        file from and directory where to save it to.
        """
        obj = str.__new__(cls)
        if url is not None and len(url) == 33:
            # looks like a Google Drive ID, turn into direct download link
            url = f"https://drive.usercontent.google.com/download?id={url}&confirm=t"
        obj._value_ = file_path.split("/")[-1]  # use file name as enum value

        obj._rel_path = file_path  # type: ignore[attr-defined] # noqa: SLF001
        obj._url = url  # type: ignore[attr-defined] # noqa: SLF001

        return obj

    def __str__(self) -> str:
        """File path associated with the file URL. Use str(DataFiles.some_key) if you
        want the absolute file path without auto-downloading the file if it doesn't
        exist yet, e.g. for use in script that generates the file in the first place.
        """
        return f"{self.base_dir}/{self._rel_path}"  # type: ignore[attr-defined]

    def __repr__(self) -> str:
        """Return enum attribute's string representation."""
        return f"{type(self).__name__}.{self.name}"

    @property
    def path(self) -> str:
        """Return the file path associated with the file URL if it exists, otherwise
        download the file first, then return the path.
        """
        key, url, rel_path = self.name, self._url, self._rel_path  # type: ignore[attr-defined]
        abs_path = f"{self.base_dir}/{rel_path}"  # type: ignore[attr-defined]
        if not os.path.isfile(abs_path):
            is_ipython = hasattr(__builtins__, "__IPYTHON__")
            # default to 'y' if not in interactive session, and user can't answer
            answer = (
                input(
                    f"{abs_path!r} associated with {key=} does not exist. Download it "
                    "now? This will cache the file for future use. [y/n] "
                )
                if is_ipython or sys.stdin.isatty()
                else "y"
            )
            if answer.lower().strip() == "y":
                if not is_ipython:
                    print(f"Downloading {key!r} from {url} to {abs_path} for caching")
                download_file(abs_path, url)
        return abs_path

    @property
    def url(self) -> str:
        """Return the URL associated with the file URL."""
        return self._url  # type: ignore[attr-defined]

    @property
    def rel_path(self) -> str:
        """Return the relative path of the file associated with the file URL."""
        return self._rel_path  # type: ignore[attr-defined]

    @property
    def label(self) -> str:
        """Return the label associated with the file URL."""
        return self.label_map.get(self.name, self.name)  # type: ignore[attr-defined]


Files.base_dir = DEFAULT_CACHE_DIR  # type: ignore[attr-defined]
Files.label_map = None  # type: ignore[attr-defined]


class DataFiles(Files):
    """Enum of data files with associated file directories and URLs."""

    mp_computed_structure_entries = (
        "mp/2023-02-07-mp-computed-structure-entries.json.gz",
        "https://figshare.com/ndownloader/files/40344436",
    )
    mp_elemental_ref_entries = (
        "mp/2023-02-07-mp-elemental-reference-entries.json.gz",
        "https://figshare.com/ndownloader/files/40387775",
    )
    mp_energies = (
        "mp/2023-01-10-mp-energies.csv.gz",
        "https://figshare.com/ndownloader/files/41296875",
    )
    mp_patched_phase_diagram = (
        "mp/2023-02-07-ppd-mp.pkl.gz",
        "https://figshare.com/ndownloader/files/40344451",
    )
    mp_trj_extxyz = (
        "mp/2023-11-22-mp-trj-extxyz-by-yuan.zip",
        "https://figshare.com/ndownloader/files/43302033",
    )
    # snapshot of every task (calculation) in MP as of 2023-03-16 (14 GB)
    all_mp_tasks = (
        "mp/2023-03-16-all-mp-tasks.zip",
        "https://figshare.com/ndownloader/files/43350447",
    )

    wbm_computed_structure_entries = (
        "wbm/2022-10-19-wbm-computed-structure-entries.json.bz2",
        "https://figshare.com/ndownloader/files/40344463",
    )
    wbm_initial_structures = (
        "wbm/2022-10-19-wbm-init-structs.json.bz2",
        "https://figshare.com/ndownloader/files/40344466",
    )
    wbm_cses_plus_init_structs = (
        "wbm/2022-10-19-wbm-computed-structure-entries+init-structs.json.bz2",
        "https://figshare.com/ndownloader/files/40344469",
    )
    wbm_summary = (
        "wbm/2023-12-13-wbm-summary.csv.gz",
        "https://figshare.com/ndownloader/files/44225498",
    )
    alignn_checkpoint = (
        "2023-06-02-pbenner-best-alignn-model.pth.zip",
        "https://figshare.com/ndownloader/files/41233560",
    )
    mp_trj = (
        "mp/2022-09-16-mp-trj-2022-09.json.gz",
        "https://figshare.com/ndownloader/files/41619375",
    )


df_wbm = pd.read_csv(DataFiles.wbm_summary.path)
# str() around Key.mat_id added for https://github.com/janosh/matbench-discovery/issues/81
df_wbm.index = df_wbm[str(Key.mat_id)]
