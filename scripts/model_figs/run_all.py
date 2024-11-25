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
scripts -= set(glob(f"{module_dir}/single_model_*.py"))  # ignore single-model scripts


# %%
exceptions: dict[str, Exception] = {}  # Collect exceptions here

for show_non_compliant in (False, True):
    for file in (pbar := tqdm(scripts)):
        pbar.set_postfix_str(file)
        init_globals = {"show_non_compliant": show_non_compliant}
        try:
            if file.endswith("tiles_energy_parity.py"):
                for which_energy in ("e-form", "each"):
                    runpy.run_path(
                        file,
                        init_globals=init_globals | {"which_energy": which_energy},
                    )
            elif file.endswith("cumulative_metrics.py"):
                for metrics in (("MAE",), ("Precision", "Recall")):
                    runpy.run_path(
                        file, init_globals=init_globals | {"metrics": metrics}
                    )
            else:
                runpy.run_path(file, init_globals=init_globals)
        except Exception as exc:
            exceptions[file] = exc

# Raise a combined exception if any errors were collected
if exceptions:
    error_messages = "\n".join(
        f"{file!r} failed: {exc!r}" for file, exc in exceptions.items()
    )
    raise RuntimeError(f"The following errors occurred:\n{error_messages}")
