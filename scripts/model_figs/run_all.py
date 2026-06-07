"""This script runs all python files in this directory which should contain all
scripts needed to generate the interactive and static PDF versions of each
model-comparison figure.
"""

# %%
import os
import runpy
import sys
from glob import glob

from pymatviz.enums import Key
from tqdm import tqdm

__author__ = "Janosh Riebesell"
__date__ = "2023-07-14"

# Inject --no-show to suppress Plotly figures opening in browser
if "--no-show" not in sys.argv:
    sys.argv.append("--no-show")

# Import cli after modifying sys.argv to apply the monkey-patch
from matbench_discovery.cli import cli_args  # noqa: F401

module_dir = os.path.dirname(__file__)

# subtract __file__ to prevent this file from calling itself
scripts = set(glob(f"{module_dir}/*.py")) - {__file__}
scripts -= set(glob(f"{module_dir}/single_model_*.py"))  # ignore single-model scripts


# %%
exceptions: dict[str, Exception] = {}  # Collect exceptions here

for show_non_compliant in (False, True):
    for script_path in (pbar := tqdm(scripts)):
        pbar.set_postfix_str(script_path)
        init_globals = {"show_non_compliant": show_non_compliant}
        try:
            if script_path.endswith("tiles_energy_parity.py"):
                # must match Key.e_form/Key.each which tiles_energy_parity.py
                # compares against (read there via globals().get("which_energy"))
                for which_energy in (Key.e_form, Key.each):
                    runpy.run_path(
                        script_path,
                        init_globals=init_globals | {"which_energy": which_energy},
                    )
            elif script_path.endswith("cumulative_metrics.py"):
                for metrics in (("MAE",), ("Precision", "Recall")):
                    runpy.run_path(
                        script_path, init_globals=init_globals | {"metrics": metrics}
                    )
            else:
                runpy.run_path(script_path, init_globals=init_globals)
        # broad on purpose: aggregate any failure from arbitrary sub-scripts run via
        # runpy and re-raise them together below (nothing is swallowed)
        except Exception as exc:  # noqa: BLE001
            exceptions[script_path] = exc

# Raise a combined exception if any errors were collected
if exceptions:
    error_messages = "\n".join(
        f"{file!r} failed: {exc!r}" for file, exc in exceptions.items()
    )
    raise RuntimeError(f"The following errors occurred:\n{error_messages}")
