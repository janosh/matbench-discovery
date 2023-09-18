"""Global variables used all across the matbench_discovery package."""

import os
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(__file__))  # repo root directory
SITE_FIGS = f"{ROOT}/site/src/figs"  # directory for interactive figures
SITE_MODELS = f"{ROOT}/site/src/routes/models"  # directory to write model analysis
FIGSHARE = f"{ROOT}/data/figshare"
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
