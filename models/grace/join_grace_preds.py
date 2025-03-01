import argparse  # Import the argparse module
import os
import re
from glob import glob

import pandas as pd
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatviz.enums import Key
from tqdm import tqdm

from matbench_discovery.data import as_dict_handler, df_wbm
from matbench_discovery.energy import calc_energy_from_e_refs, mp_elemental_ref_energies
from matbench_discovery.enums import DataFiles, MbdKey, Task

__author__ = "Yury Lysogorskiy"
__date__ = "2025-02-06"


energy_col = "grace_energy"
e_form_grace_col = "e_form_per_atom_grace"
struct_col = "grace_structure"


module_dir = os.path.dirname(__file__)
task_type = Task.IS2RE  # or Task.RS2RE, depending on what you processed


def process_results(path: str) -> None:
    """
    Processes relaxation results from a given path.

    Args:
        path (str): The path to the directory containing the .json.gz files from each
            job in the job array.
    """

    if match := re.match(r"^(\d{4}-\d{2}-\d{2})-(.*)-wbm.*", path):
        date, model_name = match[0], match[1]
        print(f"{date=} {model_name=}")
    else:
        raise ValueError(f"{path=} failed to match regex")

    glob_pattern = f"{path}/production-*.json.gz"
    file_paths = glob(glob_pattern)

    out_dir = file_paths[0].rsplit("/", 1)[0]

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
            print(f"Error reading {fn}: {e}")  # Print any errors during file reading.
            continue  # Continue to the next file

    tot_df = pd.concat(dfs)
    tot_df["id_tuple"] = (
        tot_df["material_id"].str.split("-").map(lambda x: (int(x[1]), int(x[2])))
    )
    tot_df = (
        tot_df.sort_values("id_tuple").reset_index(drop=True).drop(columns=["id_tuple"])
    )

    df_grace = tot_df.set_index("material_id")  # .drop(columns=[struct_col])

    # Create ComputedStructureEntry objects with GRACE energies and structures
    wbm_cse_path = DataFiles.wbm_computed_structure_entries.path
    df_cse = pd.read_json(wbm_cse_path).set_index(Key.mat_id)

    df_cse[Key.computed_structure_entry] = [
        ComputedStructureEntry.from_dict(dct)
        for dct in tqdm(df_cse[Key.computed_structure_entry], desc="Hydrate CSEs")
    ]

    # %% transfer ML energies and relaxed structures WBM CSEs since MP2020 energy
    # corrections applied below are structure-dependent (for oxides and sulfides)
    cse: ComputedStructureEntry
    for row in tqdm(
        df_grace.itertuples(), total=len(df_grace), desc="ML energies to CSEs"
    ):
        mat_id, struct_dict, grace_energy, *_ = row
        mlip_struct = Structure.from_dict(struct_dict)
        cse = df_cse.loc[mat_id, Key.computed_structure_entry]
        cse._energy = grace_energy  # noqa: SLF001 cse._energy is the uncorrected energy
        cse._structure = mlip_struct  # noqa: SLF001
        df_grace.loc[mat_id, Key.computed_structure_entry] = cse

    # Apply MP2020 energy corrections
    print("Applying MP2020 energy corrections")
    processed = MaterialsProject2020Compatibility().process_entries(
        df_grace[Key.computed_structure_entry], verbose=True, clean=True
    )
    if len(processed) != len(df_grace):
        raise ValueError(
            f"not all entries processed: {len(processed)=} {len(df_grace)=}"
        )

    df_grace[e_form_grace_col] = [
        calc_energy_from_e_refs(
            dict(
                composition=row["formula"],
                energy=row[Key.computed_structure_entry].energy,  # use corrected energy
            ),
            ref_energies=mp_elemental_ref_energies,
        )
        for _, row in tqdm(df_grace.iterrows(), total=len(df_grace))
    ]

    df_grace["e_form_per_atom_grace_uncorrected"] = [
        calc_energy_from_e_refs(
            dict(energy=energy, composition=formula),
            ref_energies=mp_elemental_ref_energies,
        )
        for energy, formula in tqdm(
            zip(df_grace[energy_col], df_grace["formula"]),
            total=len(df_grace),
        )
    ]

    # save relaxed structures
    df_grace.to_json(
        f"{out_dir}/{model_name}_{date}-wbm-IS2RE-FIRE.json.gz",
        default_handler=as_dict_handler,
    )
    df_grace = df_grace.drop(columns=[struct_col, Key.computed_structure_entry])

    df_wbm[[*df_grace]] = df_grace

    # %%
    bad_mask = abs(df_wbm[e_form_grace_col] - df_wbm[MbdKey.e_form_dft]) > 5
    n_preds = len(df_wbm[e_form_grace_col].dropna())
    print(f"{sum(bad_mask)=} is {sum(bad_mask) / len(df_wbm):.2%} of {n_preds:,}")

    df_grace = df_grace.round(4)
    df_grace.select_dtypes("number").to_csv(f"{out_dir}/{model_name}_{date}.csv.gz")
    df_grace.reset_index().to_json(
        f"{out_dir}/{model_name}_{date}.json.gz", default_handler=as_dict_handler
    )
    df_bad = df_grace[bad_mask].copy()
    df_bad[MbdKey.e_form_dft] = df_wbm[MbdKey.e_form_dft]
    df_bad.to_csv(f"{out_dir}/{model_name}_{date}_bad.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process relaxation results.")
    parser.add_argument(
        "path", type=str, help="Path to the directory with relaxation results."
    )

    args = parser.parse_args()
    process_results(args.path)
