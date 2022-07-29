# %%
from datetime import datetime

import pandas as pd

from ml_stability import ROOT


"""
Change JSON orientation of wbm-cleaned.json.gz and WBM step IDs to match the dielectric
Pareto frontier project.
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-06-26"

today = f"{datetime.now():%Y-%m-%d}"


def increment_wbm_material_id(wbm_id: str) -> str:
    """Turns 'wbm_step_1_0' into 'wbm_step_1_1' etc."""
    *_, step_num, material_num = wbm_id.split("_")

    assert step_num.isdigit() and material_num.isdigit()

    return f"wbm-step-{step_num}-{int(material_num) + 1}"


# %%
df_wbm = pd.read_json(f"{ROOT}/data/wbm-cleaned.json.gz", orient="split")
df_wbm.index = df_wbm.index.map(increment_wbm_material_id)

df_wbm.reset_index().to_json(
    f"{ROOT}/data/{today}-wbm-cses-and-initial-structures.json.gz"
)


# %% 2022-07-18 also increment material_ids in wbm-e-above-mp-hull.csv
for filename in (
    # "wbm-e-above-mp-hull",
    # "wren-mp-initial-structures",
    # "cgcnn-mp-initial-structures",
    # "voronoi-mp-initial-structures",
    "wren-mp-cse",
    "cgcnn-mp-cse",
    "voronoi-mp-cse",
):
    file_path = f"{ROOT}/data/2022-06-11-from-rhys/{filename}.csv"

    df = pd.read_csv(file_path)

    df["material_id"] = df.material_id.map(increment_wbm_material_id)

    df.to_csv(file_path, index=False)
