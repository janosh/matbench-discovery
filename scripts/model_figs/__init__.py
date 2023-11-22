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

__author__ = "Janosh Riebesell"
__date__ = "2023-07-14"

module_dir = os.path.dirname(__file__)

# monkey patch go.Figure.show() and Dash.run() to prevent them opening a browser
go.Figure.show = lambda self, *args, **kwargs: None
Dash.run = lambda self, *args, **kwargs: None


# %%
for file in glob(f"{module_dir}/*.py"):
    if file == __file__:  # skip this file
        continue
    print(f"Running {file.split(os.path.sep)[-1]}...")
    try:
        runpy.run_path(file)
    except Exception as exc:
        print(f"{file!r} failed: {exc}")
