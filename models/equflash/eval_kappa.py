"""Evaluate thermal conductivity predictions against DFT reference."""

# /// script
# requires-python = ">=3.11,<3.13"
# dependencies = [
# "matbench-discovery"==1.3.1,
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../../", editable = true }
# ///

import argparse

import pandas as pd
from k_srme import DFT_NONAC_REF, ID, glob2df
from k_srme.benchmark import get_metrics, process_benchmark_descriptors

parser = argparse.ArgumentParser()
parser.add_argument("--outdir", required=True, help="Output directory")
args = parser.parse_args()
model_name = "EquFlash"

json_file = f"{args.outdir}/k_srme.json.gz"
txt_path = f"{args.outdir}/metrics.txt"
in_file = f"{args.outdir}/conductivity-*.json.gz"

dft_results_file = DFT_NONAC_REF

df_mlp_results = glob2df(in_file, max_files=None).set_index(ID)
df_dft_results = pd.read_json(dft_results_file).set_index(ID)

df_mlp_filtered = df_mlp_results[df_mlp_results.index.isin(df_dft_results.index)]
df_mlp_filtered = df_mlp_filtered.reindex(df_dft_results.index)
df_mlp_processed = process_benchmark_descriptors(df_mlp_filtered, df_dft_results)

msre, msrme, _, _ = get_metrics(df_mlp_filtered)

df_mlp_processed = df_mlp_processed.round(5)
df_mlp_processed.index.name = ID
df_mlp_processed.reset_index().to_json(json_file)

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
df_mlp_print = df_mlp_filtered[["SRME", "SRE", "kappa_TOT_ave", "DFT_kappa_TOT_ave"]]

for col in ["DFT_kappa_TOT_ave", "kappa_TOT_ave"]:
    df_mlp_print[col] = df_mlp_print[col].apply(
        lambda val: val[0] if not pd.isna(val) else val
    )

df_mlp_print["SRME_failed"] = df_mlp_print["SRME"] == 2

with open(txt_path, "w") as file:
    print(f"MODEL: {model_name}", file=file)
    print(f"\tmean SRME: {msrme}", file=file)
    print(f"\tmean SRE: {msre}", file=file)
    print(df_mlp_print.round(4), file=file)
    print("\n", file=file)

    failed_ids = df_mlp_print[df_mlp_print["kappa_TOT_ave"].isna()].index.tolist()
    for mat_id in failed_ids:
        row = df_mlp_results.loc[mat_id]
        print(f"{mat_id}:\t", file=file, end="")

        error_count = 0
        if row.get("errors") and len(row["errors"]) > 0:
            print(f"errors: {row['errors']}", file=file)
            error_count += 1
        if row.get("error_traceback") and len(row["error_traceback"]) > 0:
            print(f"\t\t: {row['error_traceback']}", file=file)
            error_count += 1
        if row.get("broken_symmetry"):
            print(f"broken symmetry: {row['broken_symmetry']}", file=file)
            error_count += 1
        if row.get("imaginary_freqs"):
            print("imaginary frequencies detected", file=file)
            error_count += 1
        if error_count == 0:
            print("unknown error", file=file)

print(f"MODEL: {model_name}")
print(f"\tmean SRME: {msrme}")
print(f"\tmean SRE: {msre}")
