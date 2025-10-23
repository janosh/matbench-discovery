import argparse
import os
from glob import glob

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    args = parser.parse_args()
    files = sorted(glob(f"{args.results}/0*_*.json.gz"))
    dfs = {}
    for file_path in files:
        if file_path in dfs:
            continue
        df_i = pd.read_json(file_path).set_index("material_id")
        dfs[file_path] = df_i
    df_cat = pd.concat(dfs.values()).round(4)
    rmsd = df_cat["mlff_rmsd"].fillna(1.0).mean()

    df_7net = pd.concat(dfs.values())
    df_metadatas = pd.read_csv(
        "/home/gpu1/zetta/aixsim/1_umlff/datasets/matbench-discovery/1.0.0/wbm/wbm_formation_correction_natom.csv.gz"
    )
    df_metadatas = df_metadatas.set_index("material_id")
    df_merge = pd.concat([df_metadatas, df_7net], axis=1)
    df_merge["mlff_e_per_form"] = (
        df_merge["mlff_energy"]
        + df_merge["correction"]
        - df_merge["formation_ref_energy"]
    ) / df_merge["num_atoms"]

    df_preds = df_merge.select_dtypes("number")

    BASEDIR = "/home/gpu1/zetta/aixsim/1_umlff/datasets/matbench-discovery/1.0.0"
    df_wbm = (
        pd.read_csv(os.path.join(BASEDIR, "wbm/2023-12-13-wbm-summary.csv.gz"))
        .set_index("material_id")
        .round(3)
    )
    df_preds = df_wbm.copy()
    df_preds["e_form_per_atom_mlff"] = df_merge["mlff_e_per_form"]
    df_uniq = df_preds[df_wbm["unique_prototype"]].copy()
    df_uniq["e_form_per_atom_mlff_above_hull"] = (
        df_uniq["e_above_hull_mp2020_corrected_ppd_mp"]
        + df_uniq["e_form_per_atom_mlff"]
        - df_uniq["e_form_per_atom_mp2020_corrected"]
    )
    model_true = df_uniq["e_form_per_atom_mlff_above_hull"] <= 0
    model_false = df_uniq["e_form_per_atom_mlff_above_hull"] > 0
    actual_true = df_uniq["e_above_hull_mp2020_corrected_ppd_mp"] <= 0
    actual_false = df_uniq["e_above_hull_mp2020_corrected_ppd_mp"] > 0
    TP = model_true & actual_true
    FN = actual_true & model_false
    FP = actual_false & model_true
    n_true_pos = TP.sum()
    n_false_pos = FP.sum()
    n_false_neg = FN.sum()
    n_total_pos = n_true_pos + n_false_neg
    precision = n_true_pos / (n_true_pos + n_false_pos)  # model's discovery rate
    recall = n_true_pos / n_total_pos
    f1 = 2 * (precision * recall) / (precision + recall)

    print(f"F1\t:\t{f1:.4f}\nRMSD\t:\t{rmsd:.4f}")


if __name__ == "__main__":
    main()
