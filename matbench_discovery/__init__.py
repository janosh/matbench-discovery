"""Global variables used all across the matbench_discovery package."""

import json
import os
import warnings
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(__file__))  # repo root directory
DATA_DIR = f"{ROOT}/data"  # directory to store raw data
SITE_FIGS = f"{ROOT}/site/src/figs"  # directory for interactive figures
SITE_MODELS = f"{ROOT}/site/src/routes/models"  # directory to write model analysis
FIGSHARE = f"{ROOT}/data/figshare"
SCRIPTS = f"{ROOT}/scripts"
PDF_FIGS = f"{ROOT}/paper/figs"  # directory for light-themed PDF figures

for directory in [SITE_FIGS, SITE_MODELS, FIGSHARE, PDF_FIGS]:
    os.makedirs(directory, exist_ok=True)

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
for lineno in (120, 221, 1043):
    warnings.filterwarnings(
        action="ignore", category=UserWarning, module="pymatgen", lineno=lineno
    )

id_col = "material_id"
init_struct_col = "initial_structure"
struct_col = "structure"
e_form_col = "formation_energy_per_atom"
formula_col = "formula"
stress_col = "stress"
stress_trace_col = "stress_trace"

# load figshare 1.0.0
with open(f"{FIGSHARE}/1.0.0.json") as file:
    FIGSHARE_URLS = json.load(file)
