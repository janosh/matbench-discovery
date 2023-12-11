"""Analyze structures and composition with largest mean error across all models.
Maybe there's some chemistry/region of materials space that all models struggle with?
Might point to deficiencies in the data or models architecture.
"""


# %%
import os
from glob import glob

import numpy as np
import pandas as pd
from matminer.featurizers.site import CrystalNNFingerprint
from matminer.featurizers.structure import SiteStatsFingerprint
from pymatgen.core import Structure
from tqdm import tqdm

from matbench_discovery import DATA_DIR, id_col, timestamp
from matbench_discovery.data import DATA_FILES
from matbench_discovery.slurm import slurm_submit

__author__ = "Janosh Riebesell"
__date__ = "2023-03-26"


# %% compute all initial and final MP/WBM structure fingerprints
data_name = "wbm"
data_path = {
    "wbm": DATA_FILES.wbm_cses_plus_init_structs,
    "mp": DATA_FILES.mp_computed_structure_entries,
}[data_name]

slurm_array_task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "0"))
slurm_array_task_count = 100

out_dir = f"{DATA_DIR}/{data_name}"
os.makedirs(out_dir, exist_ok=True)

slurm_vars = slurm_submit(
    job_name=f"{data_name}-struct-fingerprints",
    out_dir=out_dir,
    partition="icelake-himem",
    account="LEE-SL3-CPU",
    time="6:0:0",
    array=f"1-{slurm_array_task_count}",
    slurm_flags=("--mem", "30G"),
)


# %%
out_path = f"{out_dir}/site-stats-{slurm_array_task_id:>03}.json.gz"
if os.path.isfile(out_path):
    raise SystemExit(f"{out_path=} already exists, exciting early")

print(f"\nJob started running {timestamp}")
print(f"{out_path=}")


# %%
df_in = pd.read_json(data_path).set_index(id_col)
if slurm_array_task_count > 1:
    df_in = np.array_split(df_in, slurm_array_task_count)[slurm_array_task_id - 1]

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
            df_in.loc[row.Index, fp_col] = ss_fp
        except Exception as exc:
            print(f"{fp_col} for {row.Index} failed: {exc}")

df_in.filter(like="site_stats_fingerprint").reset_index().to_json(out_path)


# %%
running_as_slurm_job = os.getenv("SLURM_JOB_ID")
if running_as_slurm_job:
    print(f"Job wrote {out_path=} and finished at {timestamp}")
    raise SystemExit(0)


# %%
out_files = glob(f"{out_dir}/site-stats-*.json.gz")

found_idx = [int(name.split("-")[-1].split(".")[0]) for name in out_files]
print(f"Found {len(out_files)=:,}")
missing_files = sorted(set(range(1, slurm_array_task_count + 1)) - set(found_idx))
if missing_files:
    print(f"{len(missing_files)=}: {missing_files}")

df_out = pd.concat(pd.read_json(out_file) for out_file in tqdm(out_files))
df_out = df_out.set_index(id_col)


# %%
fp_diff_col = "site_stats_fingerprint_init_final_norm_diff"
df_out[fp_diff_col] = (
    df_out[final_fp_col].map(np.array) - df_out[init_fp_col].map(np.array)
).map(np.linalg.norm)

df_out[fp_diff_col].hist(bins=100, backend="plotly")


# %%
df_out.reset_index().to_json(f"{out_dir}/site-stats.json.gz")
