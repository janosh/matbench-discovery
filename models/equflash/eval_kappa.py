import os
from glob import glob
import pandas as pd

from k_srme import glob2df, ID, DFT_NONAC_REF
from k_srme.benchmark import get_metrics, process_benchmark_descriptors

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--name", "-n", default="SevenNet")
parser.add_argument("--outdir", required="True")
args = parser.parse_args()
model_name = args.name


json_file = f"{args.outdir}/k_srme.json.gz"
txt_path = f"{args.outdir}/metrics.txt"
in_file = f"{args.outdir}/conductivity-*.json.gz"
in_folder = f"{args.outdir}"

# Comparing with the DFT results without
# non-analytical correction term (NAC)
DFT_RESULTS_FILE = DFT_NONAC_REF


module_dir = os.path.dirname(__file__)
in_pattern = f"{in_file}"
out_path = f"{json_file}"

df_mlp_results = glob2df(in_pattern, max_files=None).set_index(ID)

# df_mlp_results.reset_index().to_json(out_path)

# Read DFT results
df_dft_results = pd.read_json(DFT_RESULTS_FILE).set_index(ID)


df_mlp_filtered = df_mlp_results[df_mlp_results.index.isin(df_dft_results.index)]
df_mlp_filtered = df_mlp_filtered.reindex(df_dft_results.index)

df_mlp_processed = process_benchmark_descriptors(df_mlp_filtered, df_dft_results)

msre, msrme, _, _ = get_metrics(df_mlp_filtered)


# Save results
df_mlp_processed.round(5)
df_mlp_processed.index.name = ID
df_mlp_processed.reset_index().to_json(out_path)


# Print
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
df_mlp_print = df_mlp_filtered[
    [
        "SRME",
        "SRE",
        "kappa_TOT_ave",
        "DFT_kappa_TOT_ave",
    ]
]
df_mlp_print["DFT_kappa_TOT_ave"] = df_mlp_print["DFT_kappa_TOT_ave"].apply(
    lambda x: x[0] if not pd.isna(x) else x
)
df_mlp_print["kappa_TOT_ave"] = df_mlp_print["kappa_TOT_ave"].apply(
    lambda x: x[0] if not pd.isna(x) else x
)
df_mlp_print["SRME_failed"] = df_mlp_print["SRME"].apply(lambda x: x == 2)


with open(txt_path, "w") as f:
    print(f"MODEL: {model_name}", file=f)
    print(f"\tmean SRME: {msrme}", file=f)
    print(f"\tmean SRE: {msre}", file=f)

    print(df_mlp_print.round(4), file=f)
    print("", file=f)
    print("", file=f)
    id_index = df_mlp_print[df_mlp_print["kappa_TOT_ave"].isna()].index.tolist()
    for mat_id in id_index:
        print(mat_id, end=":\t", file=f)
        row = df_mlp_results.loc[mat_id]
        count = 0
        if len(row["errors"]) != 0:
            print(f"errors\t: {row['errors']}", file=f)
            count = count + 1
        if len(row["error_traceback"]) != 0:
            print(f"\t\t: {row['error_traceback']}", file=f)
            count = count + 1
        if row["broken_symmetry"]:
            print(f"broken\t: {row['broken_symmetry']}", file=f)
            count = count + 1
        if row["imaginary_freqs"]:
            print("imaginary frequency happened!", file=f)
            count = count + 1
        if count == 0:
            print("something unexpected happened", file=f)


print(f"MODEL: {model_name}")
print(f"\tmean SRME: {msrme}")
print(f"\tmean SRE: {msre}")
