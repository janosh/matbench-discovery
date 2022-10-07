# %%
import gzip
import os
import pickle
import warnings
from datetime import datetime
from importlib.metadata import version

import numpy as np
import pandas as pd
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from tqdm import tqdm

from mb_discovery import ROOT
from mb_discovery.slurm import slurm_submit_python

__author__ = "Janosh Riebesell"
__date__ = "2022-10-04"

"""
To slurm submit this file, run:

python path/to/file.py slurm-submit
"""

today = f"{datetime.now():%Y-%m-%d}"
module_dir = os.path.dirname(__file__)
slurm_array_task_count = 50
slurm_mem_per_node = 25_000

slurm_submit_python(
    job_name := "wbm-e-above-hull-mp",
    out_dir := f"{module_dir}/{today}-{job_name}",
    max_time := "3:0:0",
    array=f"1-{slurm_array_task_count}",
    slurm_flags=("--mem", str(slurm_mem_per_node)),
    partition="icelake-himem",
)


# %%
data_path = f"{ROOT}/data/2022-06-26-wbm-cses-and-initial-structures.json.gz"
df_wbm = pd.read_json(data_path).set_index("material_id")

slurm_job_id = os.environ.get("SLURM_JOB_ID", "debug")
slurm_array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
out_path = f"{out_dir}/{slurm_array_task_id}.csv"


print(f"Job started running {datetime.now():%Y-%m-%d@%H-%M}")
print(f"{slurm_job_id = }")
print(f"{slurm_array_task_id = }")
print(f"{data_path = }")
print(f"{out_path = }")
print(f"{version('pymatgen') = }")


if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")

df_this_job: pd.DataFrame = np.array_split(df_wbm, slurm_array_task_count)[
    slurm_array_task_id - 1
]
wbm_computed_struct_entries = [
    ComputedStructureEntry.from_dict(x) for x in tqdm(df_this_job.cse)
]

# without filter, process_entries() prints so many warnings it can crash Jupyter kernel
warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")

wbm_computed_struct_entries = MaterialsProject2020Compatibility().process_entries(
    wbm_computed_struct_entries, verbose=True, clean=True
)


# %%
ppds: dict[str, PatchedPhaseDiagram] = {}

# import io
# import urllib.request

# # 2022-01-25-ppd-mp+wbm.pkl.gz (235 MB)
# ppd_mp_wbm_url = "https://figshare.com/files/36669624"
# zipped_file = urllib.request.urlopen(ppd_mp_wbm_url)

# ppds["mp_wbm"] = pickle.load(io.BytesIO(gzip.decompress(zipped_file.read())))


# %%
with gzip.open(f"{ROOT}/data/2022-09-18-ppd-mp.pkl.gz", "rb") as zip_file:
    ppds["mp"] = pickle.load(zip_file)


# %%
for entry in tqdm(wbm_computed_struct_entries, disable=None):

    for name, ppd in ppds.items():
        e_above_hull = ppd.get_e_above_hull(entry, on_error="ignore")
        df_this_job.at[entry.entry_id, f"e_above_hull_{name}"] = e_above_hull


df_this_job.filter(like="e_above_hull_").to_csv(out_path)


# df_hull = pd.read_csv(
#     f"{ROOT}/data/2022-06-11-from-rhys/wbm-e-above-mp-hull.csv"
# ).set_index("material_id")

# df_hull["e_above_hull_mp_wbm"] = df_wbm.e_above_hull_mp_wbm

# df_hull.plot.scatter("e_above_hull_mp_wbm", "e_above_hull_mp")
