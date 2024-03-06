"""This script runs all python files in this directory which should contain all
scripts needed to generate the interactive and static PDF versions of each
model-comparison figure.
"""

# %%
import os
import runpy
from glob import glob

import plotly.graph_objects as go
from dash import Dash
from tqdm import tqdm

__author__ = "Janosh Riebesell"
__date__ = "2023-07-14"

module_dir = os.path.dirname(__file__)

# monkey patch go.Figure.show() and Dash.run() to prevent them from opening browser
go.Figure.show = lambda *_args, **_kwargs: None
Dash.run = lambda *_args, **_kwargs: None

# subtract __file__ to prevent this file from calling itself
scripts = set(glob(f"{module_dir}/*.py")) - {__file__}


# %%
for file in (pbar := tqdm(scripts)):
    pbar.set_description(file)
    try:
        if file.endswith("parity_energy_models.py"):
            for which_energy in ("each", "e-form"):
                runpy.run_path(file, init_globals={"which_energy": which_energy})
        elif file.endswith("cumulative_metrics.py"):
            for metrics in (("MAE",), ("Precision", "Recall")):
                runpy.run_path(file, init_globals={"metrics": metrics})
        elif file.endswith("rolling_mae_vs_hull_dist_wbm_batches.py"):
            runpy.run_path(file, init_globals={"models": ("CHGNet", "MACE")})
        else:
            runpy.run_path(file)
    except Exception as exc:
        print(f"{file!r} failed: {exc}")
