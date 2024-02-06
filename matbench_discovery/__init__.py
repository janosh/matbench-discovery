"""Global variables used all across the matbench_discovery package."""

from __future__ import annotations

import json
import os
import warnings
from datetime import datetime
from importlib.metadata import Distribution

import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio

from matbench_discovery.enums import (  # noqa: F401
    Key,
    Model,
    ModelType,
    Open,
    Quantity,
    Targets,
    Task,
)

pkg_name = "matbench-discovery"
direct_url = Distribution.from_name(pkg_name).read_text("direct_url.json") or "{}"
pkg_is_editable = json.loads(direct_url).get("dir_info", {}).get("editable", False)

PKG_DIR = os.path.dirname(__file__)
# repo root directory if editable install, else the pkg directory
ROOT = os.path.dirname(PKG_DIR) if pkg_is_editable else PKG_DIR
DATA_DIR = f"{ROOT}/data"  # directory to store raw data
SITE_FIGS = f"{ROOT}/site/src/figs"  # directory for interactive figures
# directory to write model analysis for website
SITE_LIB = f"{ROOT}/site/src/lib"
SCRIPTS = f"{ROOT}/scripts"  # model and date analysis scripts
PDF_FIGS = f"{ROOT}/paper/figs"  # directory for light-themed PDF figures
FIGSHARE_DIR = f"{PKG_DIR}/figshare"

for directory in (SITE_FIGS, SITE_LIB, PDF_FIGS):
    os.makedirs(directory, exist_ok=True)

os.makedirs(MP_DIR := f"{DATA_DIR}/mp", exist_ok=True)
os.makedirs(WBM_DIR := f"{DATA_DIR}/wbm", exist_ok=True)
# JSON files with URLS to data files on figshare

# directory to store model checkpoints downloaded from wandb cloud storage
CHECKPOINT_DIR = f"{ROOT}/wandb/checkpoints"
# wandb <entity>/<project name> to record new runs to
WANDB_PATH = "janosh/matbench-discovery"

# threshold on hull distance for a material to be considered stable
STABILITY_THRESHOLD = 0

timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@")[0]

# filter pymatgen warnings that spam the logs when e.g. applying corrections to
# ComputedStructureEntries or using PatchedPhaseDiagram to get e_above_hull
# warnings are:
# > No electronegativity for Ne. Setting to NaN. This has no physical meaning
# and MaterialsProject2020Compatibility to get formation energies
# > Failed to guess oxidation states for Entry
warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")


with open(f"{FIGSHARE_DIR}/1.0.0.json") as file:
    FIGSHARE_URLS = json.load(file)

# --- start global plot settings
px.defaults.labels = Quantity.dict() | Model.dict()


global_layout = dict(
    paper_bgcolor="rgba(0,0,0,0)",
    font_size=13,
    # increase legend marker size and make background transparent
    legend=dict(itemsizing="constant", bgcolor="rgba(0, 0, 0, 0)"),
)
pio.templates["global"] = dict(layout=global_layout)
pio.templates.default = "plotly_dark+global"
px.defaults.template = "plotly_dark+global"

# https://github.com/plotly/Kaleido/issues/122#issuecomment-994906924
# when seeing MathJax "loading" message in exported PDFs,
# use pio.kaleido.scope.mathjax = None


plt.rc("font", size=14)
plt.rc("legend", fontsize=16, title_fontsize=16)
plt.rc("axes", titlesize=16, labelsize=16)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rc("figure", dpi=200, titlesize=16)
plt.rcParams["figure.constrained_layout.use"] = True
# --- end global plot settings
