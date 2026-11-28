"""Global variables used all across the matbench_discovery package."""

import hashlib
import math
import os
import tempfile
import warnings
from collections.abc import Generator, Mapping
from contextlib import contextmanager
from datetime import UTC, datetime

import numpy as np

__version__ = "1.3.1"

PKG_DIR = os.path.dirname(__file__)  # Python package directory
ROOT = os.path.dirname(PKG_DIR)  # repo root directory


def repo_relative_path(file_path: str, root: str = ROOT) -> str:
    """Return file_path relative to root as POSIX, raising if it escapes root.

    Relative inputs are resolved against root (not the cwd), so an escaping path
    like ``../outside.json`` is rejected just like an external absolute one.
    """
    root = os.path.abspath(root)
    abs_path = os.path.abspath(
        file_path if os.path.isabs(file_path) else f"{root}/{file_path}"
    )
    try:
        if os.path.commonpath([root, abs_path]) == root:
            return os.path.relpath(abs_path, root).replace(os.sep, "/")
    except ValueError:  # commonpath raises across Windows drives
        pass
    raise ValueError(f"{file_path=} must be inside repo root {root!r}")


def file_digest(path: str, algorithm: str = "sha256") -> str:
    """Stream a file into a hex digest (``md5`` is only ever used for cache checks)."""
    with open(path, mode="rb") as file:
        # usedforsecurity=False keeps md5 cache checks working on FIPS-only builds
        digest = hashlib.file_digest(
            file, lambda: hashlib.new(algorithm, usedforsecurity=False)
        )
    return digest.hexdigest()


@contextmanager
def atomic_output_path(path: str, *, overwrite: bool = True) -> Generator[str]:
    """Publish a temp sibling atomically, optionally requiring a new destination."""
    parent = os.path.dirname(path) or "."
    os.makedirs(parent, exist_ok=True)
    file_descriptor, temporary_path = tempfile.mkstemp(
        prefix=f".{os.path.basename(path)}.", suffix=".tmp", dir=parent
    )
    os.close(file_descriptor)
    try:
        yield temporary_path
        if overwrite:
            os.replace(temporary_path, path)
        else:
            os.link(temporary_path, path)
    finally:
        if os.path.isfile(temporary_path):
            os.remove(temporary_path)


def json_ready(value: object, *, nan_to_none: bool = True) -> object:
    """Recursively convert NumPy scalars and arrays, tuples and mappings to strict
    JSON data. Non-finite floats become None, or raise ValueError when nan_to_none
    is False (e.g. run events, where NaN signals a bug rather than missing data).
    """
    if isinstance(value, Mapping):
        return {
            str(key): json_ready(item, nan_to_none=nan_to_none)
            for key, item in value.items()
        }
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist(), nan_to_none=nan_to_none)  # ty: ignore[no-matching-overload]
    if isinstance(value, np.generic):
        return json_ready(value.item(), nan_to_none=nan_to_none)
    if isinstance(value, tuple | list):
        return [json_ready(item, nan_to_none=nan_to_none) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        if nan_to_none:
            return None
        raise ValueError(f"Out of range float values are not JSON compliant: {value}")
    if isinstance(value, str | int | float | bool) or value is None:
        return value
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


TEST_FILES = f"{ROOT}/tests/files"  # directory to store test data
# directory for data-only figure payloads (gzipped JSON written by
# matbench_discovery.figs, imported by Svelte pages and rendered with matterviz)
SITE_FIG_DATA = f"{ROOT}/site/src/figs"
# directory to write model analysis for website
SITE_DIR = f"{ROOT}/site/src"
PDF_FIGS = f"{ROOT}/paper/figs"  # directory for light-themed PDF figures

# directory to cache downloaded data files, the same for repo clones and pip installs
DEFAULT_CACHE_DIR = os.getenv(
    "MBD_CACHE_DIR", os.path.expanduser("~/.cache/matbench-discovery")
)

# no import-time makedirs: writers create the directories they write to
MP_DIR = f"{DEFAULT_CACHE_DIR}/mp"
WBM_DIR = f"{DEFAULT_CACHE_DIR}/wbm"

# threshold on hull distance for a material to be considered stable
STABILITY_THRESHOLD = 0

timestamp = f"{datetime.now(tz=UTC):%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@", maxsplit=1)[0]

# filter pymatgen warnings that spam the logs when e.g. applying corrections to
# ComputedStructureEntries or using PatchedPhaseDiagram to get e_above_hull
# and MaterialsProject2020Compatibility to get formation energies. Warnings are:
# > No Pauling electronegativity for Ne. Setting to NaN. This has no physical meaning
# > Failed to guess oxidation states for Entry
# pymatgen warns with stacklevel=2, attributing Element.X's warning to its caller
# (functools.cached_property, or user code), so that filter can't pin a module.
for msg, category, module in {
    ("No Pauling electronegativity for", UserWarning, ""),
    ("Failed to guess oxidation states for Entry", UserWarning, "pymatgen"),
    ("torch.load", FutureWarning, ""),
    ("logm result may be inaccurate, approximate err", RuntimeWarning, ""),
}:
    warnings.filterwarnings(
        action="ignore", category=category, module=module, message=msg
    )
