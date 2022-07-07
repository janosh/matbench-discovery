# %%
from __future__ import annotations

import warnings
from datetime import datetime
from typing import Any

import pandas as pd
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


# %%
df_wbm = pd.read_json(
    f"{ROOT}/data/2022-06-26-wbm-cses-and-initial-structures.json.gz"
).set_index("material_id")


# %%
relax_results = []
relaxer = Relaxer()  # This loads the default pre-trained M3GNet model

try:
    for row in tqdm(df_wbm.itertuples(), total=len(df_wbm)):
        init_struct = row.initial_structure
        pmg_struct = Structure.from_dict(init_struct)
        relax_result = relaxer.relax(pmg_struct)
        relax_dict = {
            "material_id": row.Index,
            "final_structure": relax_result["final_structure"],
            "trajectory": relax_result["trajectory"].__dict__,
        }
        relax_dict["trajectory"].pop("atoms")
        relax_results.append(relax_dict)
except KeyboardInterrupt:
    pass  # make long-running loop ctrl+c interruptible


# %%
def default_handler(obj: Any) -> dict[str, Any] | None:
    try:
        return obj.as_dict()
    except AttributeError:
        return None  # replace ASE atoms with None since they aren't JSON serializable


df_relax_results = pd.DataFrame(relax_results)
df_relax_results.to_json(
    f"{ROOT}/data/{timestamp}-m3gnet_wbm_relax_results.json.gz",
    default_handler=default_handler,
)


# df_results = pd.read_json(
#     f"{ROOT}/data/{timestamp}-m3gnet_wbm_relax_results.json.gz"
# ).set_index("material_id")

# pd.json_normalize(df_results.trajectory)
