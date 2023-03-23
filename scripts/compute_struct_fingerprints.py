"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""


# %%
import os
import warnings

import numpy as np
import pandas as pd
from matminer.featurizers.site import CrystalNNFingerprint
from matminer.featurizers.structure import SiteStatsFingerprint
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import ROOT, timestamp
from matbench_discovery.data import DATA_FILES
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2023-03-26"

warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")


# %% compute all initial and final MP/WBM structure fingerprints
data_name = "wbm"
data_path = {
    "wbm": DATA_FILES.wbm_cses_plus_init_structs,
    "mp": DATA_FILES.mp_computed_structure_entries,
}[data_name]

slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", 0))
slurm_array_task_count = 100

job_name = f"make-{data_name}-struct-fingerprints"
out_dir = f"{ROOT}/data/{data_name}/structure-fingerprints"
os.makedirs(out_dir, exist_ok=True)

slurm_vars = slurm_submit(
    job_name=job_name,
    out_dir=out_dir,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time="6:0:0",
    array=f"1-{slurm_array_task_count}",
)


# %%
out_path = f"{out_dir}/site-stats-{slurm_array_task_id}.json.gz"
if os.path.isfile(out_path):
    raise SystemExit(f"{out_path = } already exists, exciting early")

print(f"\nJob started running {timestamp}")
print(f"{out_path=}")


# %%
df_in: pd.DataFrame = np.array_split(
    pd.read_json(data_path).set_index("material_id"), slurm_array_task_count
)[slurm_array_task_id - 1]

cnn_fp = CrystalNNFingerprint.from_preset("ops")
# including "minimum" and "maximum" increases the fingerprint length from 61 to 122
site_stats_fp = SiteStatsFingerprint(
    cnn_fp, stats=("mean", "std_dev", "minimum", "maximum")
)


# %%
init_struct_col = "initial_structure"
final_struct_col = "computed_structure_entry"
init_fp_col = "initial_site_stats_fingerprint"
final_fp_col = "final_site_stats_fingerprint"
for struct_col, fp_col in (
    (init_struct_col, init_fp_col),
    (final_struct_col, final_fp_col),
    ("entry", final_fp_col),
):
    if struct_col not in df_in:
        continue
    df_in[fp_col] = None

    for row in tqdm(df_in.itertuples(), total=len(df_in)):
        struct = getattr(row, struct_col)
        if "structure" in struct:  # is a ComputedStructureEntry as dict
            struct = struct["structure"]
        struct = Structure.from_dict(struct)
        try:
            ss_fp = site_stats_fp.featurize(struct)
            df_in.at[row.Index, fp_col] = ss_fp
        except Exception as exc:
            print(f"{fp_col} for {row.Index} failed: {exc}")

df_in.filter(like="site_stats_fingerprint").to_json(out_path)
