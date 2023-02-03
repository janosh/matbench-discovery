"""Global variables used all across the matbench_discovery package."""

import json
import os
import sys
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(__file__))  # repo root
FIGS = f"{ROOT}/site/src/figs"  # directory to store interactive figures
STATIC = f"{ROOT}/site/static/figs"  # directory to store static figures, is symlinked
# into site/src/routes/paper/figs dir
MODELS = f"{ROOT}/site/src/routes/models"  # directory to write model analysis
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

# load docs, repo, package URLs from package.json
with open(f"{ROOT}/site/package.json") as file:
    pkg = json.load(file)
    pypi_keys_to_npm = dict(Docs="homepage", Repo="repository", Package="package")
    URLs = {key: pkg[val] for key, val in pypi_keys_to_npm.items()}
