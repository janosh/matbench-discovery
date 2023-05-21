"""Global variables used all across the matbench_discovery package."""

import os
import sys
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(__file__))  # repo root directory
FIGS = f"{ROOT}/site/src/figs"  # directory for interactive figures
MODELS = f"{ROOT}/site/src/routes/models"  # directory to write model analysis
FIGSHARE = f"{ROOT}/data/figshare"
PDF_FIGS = f"{ROOT}/paper/figs"  # directory for light-themed PDF figures

for directory in [FIGS, MODELS, FIGSHARE, PDF_FIGS]:
    os.makedirs(directory, exist_ok=True)

# whether a currently running slurm job is in debug mode
DEBUG = "DEBUG" in os.environ or (
    "slurm-submit" not in sys.argv and "SLURM_JOB_ID" not in os.environ
)
# directory to store model checkpoints downloaded from wandb cloud storage
CHECKPOINT_DIR = f"{ROOT}/wandb/checkpoints"
# wandb <entity>/<project name> to record new runs to
WANDB_PATH = "janosh/matbench-discovery"

timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@")[0]
