"""Global variables used all across the matbench_discovery package."""

import os
import warnings
from datetime import UTC, datetime

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


DATA_DIR = f"{ROOT}/data"  # directory to store raw data
TEST_FILES = f"{ROOT}/tests/files"  # directory to store test data
# directory for data-only figure payloads (gzipped JSON written by
# matbench_discovery.figs, imported by Svelte pages and rendered with matterviz)
SITE_FIG_DATA = f"{ROOT}/site/src/figs"
# directory to write model analysis for website
SITE_DIR = f"{ROOT}/site/src"
PDF_FIGS = f"{ROOT}/paper/figs"  # directory for light-themed PDF figures

# directory to cache downloaded data files
DEFAULT_CACHE_DIR = os.getenv(
    "MBD_CACHE_DIR",
    DATA_DIR  # use DATA_DIR to locally cache data files if full repo was cloned
    if os.path.isdir(DATA_DIR)
    # use ~/.cache if matbench-discovery was installed from PyPI
    else os.path.expanduser("~/.cache/matbench-discovery"),
)

for directory in (SITE_FIG_DATA, SITE_DIR, PDF_FIGS):
    os.makedirs(directory, exist_ok=True)

os.makedirs(MP_DIR := f"{DATA_DIR}/mp", exist_ok=True)
os.makedirs(WBM_DIR := f"{DATA_DIR}/wbm", exist_ok=True)
# JSON files with URLS to data files on figshare

# threshold on hull distance for a material to be considered stable
STABILITY_THRESHOLD = 0

timestamp = f"{datetime.now(tz=UTC):%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@", maxsplit=1)[0]

# filter pymatgen warnings that spam the logs when e.g. applying corrections to
# ComputedStructureEntries or using PatchedPhaseDiagram to get e_above_hull
# warnings are:
# > No electronegativity for Ne. Setting to NaN. This has no physical meaning
# and MaterialsProject2020Compatibility to get formation energies
# > Failed to guess oxidation states for Entry
for msg, category, module in {
    ("No electronegativity for", UserWarning, "pymatgen"),
    ("Failed to guess oxidation states for Entry", UserWarning, "pymatgen"),
    ("torch.load", FutureWarning, ""),
    ("logm result may be inaccurate, approximate err", RuntimeWarning, ""),
}:
    warnings.filterwarnings(
        action="ignore", category=category, module=module, message=msg
    )
