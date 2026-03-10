"""Calculate discovery metrics (F1, precision, recall) for model predictions."""

# /// script
# requires-python = ">=3.11,<3.13"
# dependencies = [
# "matbench-discovery==1.3.1",
# ]
#
# [tool.uv.sources]
# matbench-discovery = { path = "../../", editable = true }
# ///

import argparse
import os
from glob import glob

import pandas as pd


def main() -> None:
    """Calculate F1 score and RMSD for matbench-discovery predictions."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, help="Results directory")
    parser.add_argument(
        "--dataset-root",
        default=os.environ.get(
            "DISCOVERY_DATASET_ROOT", "data/matbench-discovery/1.0.0"
        ),
        help="Root dir containing wbm CSVs",
    )
    parser.add_argument(
        "--wbm-summary",
        default="wbm/2023-12-13-wbm-summary.csv.gz",
        help="Path relative to --dataset-root for WBM summary CSV",
    )
    args = parser.parse_args()
    files = sorted(glob(f"{args.results}/0*_*.json.gz"))
    dataframes = [
        pd.read_json(file_path).set_index("material_id") for file_path in files
    ]
    df_model = pd.concat(dataframes).round(4)
    rmsd = df_model["mlff_rmsd"].fillna(1.0).mean()

    meta_csv = f"{args.dataset_root}/wbm/wbm_formation_correction_natom.csv.gz"
    df_metadata = pd.read_csv(meta_csv).set_index("material_id")
    df_merge = pd.concat([df_metadata, df_model], axis=1)
    df_merge["mlff_e_per_form"] = (
        df_merge["mlff_energy"]
        + df_merge["correction"]
        - df_merge["formation_ref_energy"]
    ) / df_merge["num_atoms"]

    df_wbm = pd.read_csv(os.path.join(args.dataset_root, args.wbm_summary)).set_index(
        "material_id"
    )

    df_preds = df_wbm.copy()
    df_preds["e_form_per_atom_mlff"] = df_merge["mlff_e_per_form"]
    df_uniq = df_preds[df_wbm["unique_prototype"]].copy()
    df_uniq["e_form_per_atom_mlff_above_hull"] = (
        df_uniq["e_above_hull_mp2020_corrected_ppd_mp"]
        + df_uniq["e_form_per_atom_mlff"]
        - df_uniq["e_form_per_atom_mp2020_corrected"]
    )

    model_pred_stable = df_uniq["e_form_per_atom_mlff_above_hull"] <= 0
    actual_stable = df_uniq["e_above_hull_mp2020_corrected_ppd_mp"] <= 0

    n_true_pos = (model_pred_stable & actual_stable).sum()
    n_false_pos = (model_pred_stable & ~actual_stable).sum()
    n_false_neg = (~model_pred_stable & actual_stable).sum()

    precision = (
        n_true_pos / (n_true_pos + n_false_pos)
        if (n_true_pos + n_false_pos) > 0
        else 0.0
    )
    recall = (
        n_true_pos / (n_true_pos + n_false_neg)
        if (n_true_pos + n_false_neg) > 0
        else 0.0
    )
    f1_score = (
        2 * (precision * recall) / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )

    print(f"F1\t:\t{f1_score:.4f}\nRMSD\t:\t{rmsd:.4f}")


if __name__ == "__main__":
    main()
