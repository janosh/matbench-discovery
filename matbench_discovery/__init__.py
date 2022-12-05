from __future__ import annotations

import os
import sys
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(__file__))
DEBUG = "slurm-submit" not in sys.argv and "SLURM_JOB_ID" not in os.environ
CHECKPOINT_DIR = f"{ROOT}/wandb/checkpoints"

timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@")[0]
