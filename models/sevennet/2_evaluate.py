import os
from glob import glob

import pandas as pd
from k_srme import DFT_NONAC_REF, ID, glob2df
from k_srme.benchmark import get_metrics, process_benchmark_descriptors

model_name = "SevenNet_l3i5"

json_file = "k_srme.json.gz"
txt_path = "metrics.txt"
in_file = "conductivity_*-???.json.gz"

in_folder = f"2024-12-03-{model_name}-phononDB-LTC-FIRE_2SR_force0.0001_sym1e-05/"  # date update

# Comparing with the DFT results without
# non-analytical correction term (NAC)
DFT_RESULTS_FILE = DFT_NONAC_REF


module_dir = os.path.dirname(__file__)
in_pattern = f"{module_dir}/{in_folder}/{in_file}"
out_path = f"{module_dir}/{in_folder}/{json_file}"


# Read MLP results
if not glob(in_pattern):
    if os.path.exists(out_path):
        df_mlp_results = pd.read_json(out_path).set_index(ID)
else:
    df_mlp_results = glob2df(in_pattern, max_files=None).set_index(ID)


# df_mlp_results.reset_index().to_json(out_path)

# Read DFT results
df_dft_results = pd.read_json(DFT_RESULTS_FILE).set_index(ID)


df_mlp_filtered = df_mlp_results[df_mlp_results.index.isin(df_dft_results.index)]
df_mlp_filtered = df_mlp_filtered.reindex(df_dft_results.index)

df_mlp_processed = process_benchmark_descriptors(df_mlp_filtered, df_dft_results)

mSRE, mSRME, rmseSRE, rmseSRME = get_metrics(df_mlp_filtered)


# Save results
df_mlp_processed.round(5)
df_mlp_processed.index.name = ID
df_mlp_processed.reset_index().to_json(out_path)


# Print
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
df_mlp_print = df_mlp_filtered[
    [
        "desc",
        "SRME",
        "SRE",
        "kappa_TOT_ave",
        "DFT_kappa_TOT_ave",
        "imaginary_freqs",
        "errors",
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
    print(f"\tmean SRME: {mSRME}", file=f)
    print(f"\tmean SRE: {mSRE}", file=f)

    print(df_mlp_print.round(4), file=f)


df_mlp_print = df_mlp_print[
    ["desc", "SRME", "SRE", "kappa_TOT_ave", "DFT_kappa_TOT_ave"]
]
print(df_mlp_print.round(3))

print(f"MODEL: {model_name}")
print(f"\tmean SRME: {mSRME}")
print(f"\tmean SRE: {mSRE}")
