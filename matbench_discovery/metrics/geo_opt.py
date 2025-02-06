"""Functions to calculate and save geometry optimization metrics."""

import pandas as pd
from pymatviz.enums import Key, Task
from ruamel.yaml.comments import CommentedMap

from matbench_discovery.data import Model, round_trip_yaml
from matbench_discovery.enums import MbdKey


def write_geo_opt_metrics_to_yaml(
    df_geo_opt: pd.DataFrame, model: Model, symprec: float
) -> None:
    """Write geometry optimization metrics to model YAML metadata files.

    Args:
        df_geo_opt (pd.DataFrame): DataFrame with geometry optimization metrics as
            columns, including:
            - structure_rmsd_vs_dft: RMSD between predicted and DFT structures
            - n_sym_ops_mae: Mean absolute error in number of symmetry operations
            - symmetry_decrease: Fraction of structures with decreased symmetry
            - symmetry_match: Fraction of structures with matching symmetry
            - symmetry_increase: Fraction of structures with increased symmetry
            - n_structs: Number of structures evaluated
        model (Model): Instance of Model enum that was analyzed in df_geo_opt.
        symprec (float): symmetry precision for comparing ML and DFT relaxed structures.
    """
    # Load existing metadata
    with open(model.yaml_path) as file:
        model_metadata = round_trip_yaml.load(file)

    all_metrics = model_metadata.setdefault("metrics", {})

    # Get metrics for this model
    new_metrics = {
        str(Key.rmsd): float(round(df_geo_opt[MbdKey.structure_rmsd_vs_dft], 4)),
        str(Key.n_sym_ops_mae): float(round(df_geo_opt[Key.n_sym_ops_mae], 4)),
        str(Key.symmetry_decrease): float(round(df_geo_opt[Key.symmetry_decrease], 4)),
        str(Key.symmetry_match): float(round(df_geo_opt[Key.symmetry_match], 4)),
        str(Key.symmetry_increase): float(round(df_geo_opt[Key.symmetry_increase], 4)),
        str(Key.n_structures): int(df_geo_opt[Key.n_structures]),
    }
    symprec_key = f"{symprec=:.0e}".replace("e-0", "e-")

    geo_opt_metrics = CommentedMap(all_metrics.setdefault(Task.geo_opt, {}))
    metrics_for_symprec = CommentedMap(geo_opt_metrics.setdefault(symprec_key, {}))
    metrics_for_symprec.update(new_metrics)

    # Define units for metrics
    metric_units = {
        Key.rmsd: "Å",
        Key.n_sym_ops_mae: "unitless",
        Key.symmetry_decrease: "fraction",
        Key.symmetry_match: "fraction",
        Key.symmetry_increase: "fraction",
        Key.n_structures: "count",
    }

    # Add units as YAML end-of-line comments
    for key in new_metrics:
        if unit := metric_units.get(key):
            metrics_for_symprec.yaml_add_eol_comment(unit, key, column=1)

    geo_opt_metrics[symprec_key] = metrics_for_symprec
    all_metrics[Task.geo_opt] = geo_opt_metrics

    # Write back to file
    with open(model.yaml_path, mode="w") as file:
        round_trip_yaml.dump(model_metadata, file)


def calc_geo_opt_metrics(df_model_analysis: pd.DataFrame) -> pd.DataFrame:
    """Calculate geometry optimization metrics for a single model.

    Args:
        df_model_analysis (pd.DataFrame): DataFrame with geometry optimization metrics
            for one model. Required columns are:
            - structure_rmsd_vs_dft: RMSD between predicted and DFT structures
            - n_sym_ops_diff: Difference in number of symmetry operations vs DFT
            - spg_num_diff: Difference in space group number vs DFT
        model_name (str): Name of the model being analyzed.

    Returns:
        pd.DataFrame: DataFrame with geometry optimization metrics.
        Shape = (1, n_metrics). Columns include:
        - structure_rmsd_vs_dft: Mean RMSD between predicted and DFT structures
        - n_sym_ops_mae: Mean absolute error in number of symmetry operations
        - symmetry_decrease: Fraction of structures with decreased symmetry
        - symmetry_match: Fraction of structures with matching symmetry
        - symmetry_increase: Fraction of structures with increased symmetry
        - n_structs: Number of structures evaluated
    """
    # Get relevant columns
    spg_diff = df_model_analysis[MbdKey.spg_num_diff]
    n_sym_ops_diff = df_model_analysis[MbdKey.n_sym_ops_diff]
    rmsd = df_model_analysis[MbdKey.structure_rmsd_vs_dft]

    # Count total number of structures (excluding NaN values)
    total = len(spg_diff.dropna())

    # Calculate RMSD and MAE metrics
    mean_rmsd = rmsd.mean()
    sym_ops_mae = n_sym_ops_diff.abs().mean()

    # Count cases where spacegroup changed
    changed_mask = spg_diff != 0
    # Among changed cases, count whether symmetry increased or decreased
    sym_decreased = (n_sym_ops_diff < 0) & changed_mask
    sym_increased = (n_sym_ops_diff > 0) & changed_mask
    sym_matched = ~changed_mask

    return {
        str(MbdKey.structure_rmsd_vs_dft): float(mean_rmsd),
        str(Key.n_sym_ops_mae): float(sym_ops_mae),
        str(Key.symmetry_decrease): float(sym_decreased.sum() / total),
        str(Key.symmetry_match): float(sym_matched.sum() / total),
        str(Key.symmetry_increase): float(sym_increased.sum() / total),
        str(Key.n_structures): total,
    }
