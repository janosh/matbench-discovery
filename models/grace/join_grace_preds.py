import glob
import os
import argparse  # Import the argparse module
import re

import pandas as pd
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import as_dict_handler, df_wbm
from matbench_discovery.energy import calc_energy_from_e_refs, mp_elemental_ref_energies
from matbench_discovery.enums import MbdKey, Task

__author__ = "Yury Lysogorskiy"
__date__ = "2025-02-06"


energy_column = "grace_energy"  # or your actual column name.
e_form_grace_col = "e_form_per_atom_grace"  # or your desired name
struct_col = "grace_structure"


module_dir = os.path.dirname(__file__)
task_type = Task.IS2RE  # or Task.RS2RE, depending on what you processed


def split_string(input_string):
    """
    Splits the input string into date and model name.

    Args:
        input_string: The string to split.

    Returns:
        A tuple containing the date string and the model name string,
        or (None, None) if the string doesn't match the expected pattern.
    """
    match = re.match(r"^(\d{4}-\d{2}-\d{2})-(.*)-wbm.*", input_string)
    if match:
        date_str = match.group(1)
        model_name = match.group(2)
        return date_str, model_name
    else:
        return None, None

def process_results(path: str):
    """
    Processes relaxation results from a given path.

    Args:
        path (str): The path to the directory containing the .json.gz files.
    """
    date, model_name = split_string(path)
    print("date:",date)
    print("model_name:", model_name)

    glob_pattern = os.path.join(path, "production-*.json.gz")
    file_paths = glob.glob(glob_pattern)

    print(f"Found {len(file_paths):,} files for {glob_pattern = }")

    if not file_paths:
        print(f"No files found matching {glob_pattern}. Exiting.")
        return  # Exit if no files are found

    dfs: list[pd.DataFrame] = []
    for fn in file_paths:
        print(fn)
        try:
            dfs.append(pd.read_json(fn))
        except Exception as e:
            print(f"Error reading {fn}: {e}") # Print any errors during file reading.
            continue # Continue to the next file

    tot_df = pd.concat(dfs)
    tot_df["id_tuple"] = (
        tot_df["material_id"].str.split("-").map(lambda x: (int(x[1]), int(x[2])))
    )
    tot_df = (
        tot_df.sort_values("id_tuple")
        .reset_index(drop=True)
        .drop(columns=["id_tuple", struct_col])
    )

    df_grace = tot_df.set_index("material_id")
    df_grace[Key.formula] = df_wbm[Key.formula]


    print("Calculating formation energies")
    e_form_list = []
    for _, row in tqdm(df_grace.iterrows(), total=len(df_grace)):
        e_form = calc_energy_from_e_refs(
            row["formula"],
            ref_energies=mp_elemental_ref_energies,
            total_energy=row[energy_column],
        )
        e_form_list.append(e_form)


    df_grace[e_form_grace_col] = e_form_list

    df_wbm[[*df_grace]] = df_grace


    # %%
    bad_mask = abs(df_wbm[e_form_grace_col] - df_wbm[MbdKey.e_form_dft]) > 5
    n_preds = len(df_wbm[e_form_grace_col].dropna())
    print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")
    out_path = file_paths[0].rsplit("/", 1)[0]  # Get directory from first file path.

    df_grace = df_grace.round(4)
    df_grace.select_dtypes("number").to_csv(f"{out_path}/{model_name}_{date}.csv.gz") #added model and date
    df_grace.reset_index().to_json(f"{out_path}/{model_name}_{date}.json.gz", default_handler=as_dict_handler) #added model and date
    df_bad = df_grace[bad_mask].copy()
    df_bad[MbdKey.e_form_dft] = df_wbm[MbdKey.e_form_dft]
    df_bad.to_csv(f"{out_path}/{model_name}_{date}_bad.csv") #added model and date



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process relaxation results.")
    parser.add_argument("path", type=str, help="Path to the directory with relaxation results.")

    args = parser.parse_args()
    process_results(args.path)
