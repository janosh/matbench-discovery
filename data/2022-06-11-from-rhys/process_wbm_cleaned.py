from datetime import datetime

import pandas as pd

from ml_stability import ROOT


"""
Change JSON orientation of wbm-cleaned.json.gz and WBM step IDs to match the dielectric
Pareto frontier project.
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-06-26"


df_wbm = pd.read_json(f"{ROOT}/data/wbm-cleaned.json.gz", orient="split")


def increment_last_idx(idx: str) -> str:
    """Turns 'step-1-0' into 'wbm-step-1-1' etc."""
    parts = idx.split("_")
    incremented = int(parts[-1]) + 1
    parts[-1] = str(incremented)
    new_id = "wbm_" + "_".join(parts)
    return new_id


df_wbm.index = df_wbm.index.map(increment_last_idx)

today = f"{datetime.now():%Y-%m-%d}"
df_wbm.reset_index().to_json(
    f"{ROOT}/data/{today}-wbm-cses-and-initial-structures.json.gz"
)
