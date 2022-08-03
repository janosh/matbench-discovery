# %%
from __future__ import annotations

import warnings
from datetime import datetime
from typing import Any

import m3gnet
import pandas as pd
import tensorflow as tf
from m3gnet.models import Relaxer
from pymatgen.core import Structure
from tqdm import tqdm

from ml_stability import ROOT


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


df_relax_results = pd.DataFrame(relax_results).T
df_relax_results.index.name = "material_id"
df_relax_results.to_json(
    f"{ROOT}/data/{timestamp}-m3gnet-wbm-relax-results.json.gz",
    default_handler=default_handler,
)


# %% merge previous and new results
df_results_old = pd.read_json(
    f"{ROOT}/data/2022-07-17@20-27-m3gnet-wbm-relax-results.json.gz"
).set_index("material_id")
df_results_new = pd.read_json(
    f"{ROOT}/data/2022-08-02@16-59-m3gnet-wbm-relax-results.json.gz"
).set_index("material_id")

df_results_new = pd.concat([df_results_old, df_results_new])

df_results_new.index.name = "material_id"


# WARNING: overwrites original file, make sure new df is as desired
df_results_new.reset_index().to_json(
    f"{ROOT}/data/2022-08-02@16-59-m3gnet-wbm-relax-results.json.gz"
)
