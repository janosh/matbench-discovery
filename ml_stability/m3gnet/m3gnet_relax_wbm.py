# %%
from __future__ import annotations

import warnings
from datetime import datetime
from typing import Any

import m3gnet
import pandas as pd
import tensorflow as tf
from diel_frontier.patched_phase_diagram import load_ppd
from m3gnet.models import Relaxer
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core import Structure
from tqdm import tqdm

from ml_stability import ROOT
from ml_stability.plots.plot_funcs import hist_classify_stable_as_func_of_hull_dist


"""
To slurm submit this file, prepend its invocation with

```sh
. /etc/profile.d/modules.sh
module load rhel8/default-amp cuda/11.2 cudnn
```

Takes about 40h for full WBM dataset on 1 A100 GPU.
1%: 1360/257486 [20:06<32:57:32,  2.16it/s]

Requires regular MEGNet pip installation:
pip install m3gnet
"""

__author__ = "Janosh Riebesell"
__date__ = "2022-06-18"

warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")
warnings.filterwarnings(action="ignore", category=UserWarning, module="tensorflow")

timestamp = f"{datetime.now():%Y-%m-%d@%H-%M}"
relax_results: dict[str, dict[str, Any]] = {}

print(f"Using M3GNet {m3gnet.__version__} version")


# %%
df_wbm = pd.read_json(
    f"{ROOT}/data/2022-06-26-wbm-cses-and-initial-structures.json.gz"
).set_index("material_id")


# %%
tf.config.list_physical_devices()


# %%
relaxer = Relaxer()  # This loads the default pre-trained M3GNet model

try:
    for material_id, init_struct in tqdm(
        df_wbm.initial_structure.items(), total=len(df_wbm)
    ):
        if material_id in relax_results:
            continue
        pmg_struct = Structure.from_dict(init_struct)
        relax_result = relaxer.relax(pmg_struct)
        relax_dict = {
            "final_structure": relax_result["final_structure"],
            "trajectory": relax_result["trajectory"].__dict__,
        }
        # remove non-serializable AseAtoms from trajectory
        relax_dict["trajectory"].pop("atoms")
        relax_results[material_id] = relax_dict

except KeyboardInterrupt:
    print("Interrupted")  # make long-running loop ctrl+c interruptible


# %%
def default_handler(obj: Any) -> dict[str, Any] | None:
    try:
        return obj.as_dict()
    except AttributeError:
        return None  # replace ASE atoms with None since they aren't JSON serializable


df_m3gnet = pd.DataFrame(relax_results).T
df_m3gnet.index.name = "material_id"
df_m3gnet.to_json(
    f"{ROOT}/data/{timestamp}-m3gnet-wbm-relax-results.json.gz",
    default_handler=default_handler,
)


# %% 2022-08-03 --- merge previous and new results
df_results_old = pd.read_json(
    f"{ROOT}/data/2022-07-17@20-27-m3gnet-wbm-relax-results.json.gz"
).set_index("material_id")
df_m3gnet = pd.read_json(
    f"{ROOT}/data/2022-08-02@16-59-m3gnet-wbm-relax-results.json.gz"
).set_index("material_id")

df_m3gnet = pd.concat([df_results_old, df_m3gnet])

df_m3gnet.index.name = "material_id"


# WARNING: overwrites original file, make sure new df is as desired
out_file = "2022-08-02@16-59-m3gnet-wbm-relax-results-with-e_form-and-pd_entry.json.gz"
df_m3gnet.reset_index().to_json(f"{ROOT}/data/{out_file}")


# %%
df_m3gnet["m3gnet_structure"] = df_m3gnet.final_structure.map(Structure.from_dict)
df_m3gnet["m3gnet_energy"] = df_m3gnet.trajectory.map(lambda x: x["energies"][-1][0])


ppd_mp_wbm = load_ppd("ppd-mp+wbm-2022-01-25.pkl.gz")


df_m3gnet["pd_entry"] = [
    PDEntry(row.m3gnet_structure.composition, row.m3gnet_energy)
    for row in df_m3gnet.itertuples()
]
df_m3gnet["e_form_m3gnet"] = df_m3gnet.pd_entry.map(ppd_mp_wbm.get_form_energy_per_atom)


df_m3gnet.hist(bins=80, figsize=(22, 5), layout=(1, 3))
df_m3gnet.isna().sum()


# %%
df_hull = pd.read_csv(
    f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
).set_index("material_id")

df_m3gnet["e_above_mp_hull"] = df_hull.e_above_mp_hull

df_summary = pd.read_csv(f"{ROOT}/data/wbm-steps-summary.csv", comment="#").set_index(
    "material_id"
)

df_m3gnet["e_form_wbm"] = df_summary.e_form


# %%
df_m3gnet.hist(bins=80, figsize=(18, 12))
df_m3gnet.isna().sum()


# %%
ax_hull_dist_hist = hist_classify_stable_as_func_of_hull_dist(
    formation_energy_targets=df_m3gnet.e_form_wbm,
    formation_energy_preds=df_m3gnet.e_form_m3gnet,
    e_above_hull_vals=df_m3gnet.e_above_mp_hull,
)

ax_hull_dist_hist.figure.savefig(
    f"{ROOT}/data/2022-08-02@16-59-m3gnet-wbm-hull-dist-hist.pdf"
)
