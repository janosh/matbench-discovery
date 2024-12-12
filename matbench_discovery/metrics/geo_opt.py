"""Functions to calculate and save geometry optimization metrics."""

import pandas as pd
from pymatviz.enums import Key, Task
from ruamel.yaml.comments import CommentedMap

from matbench_discovery.data import Model, round_trip_yaml
from matbench_discovery.enums import MbdKey


def write_geo_opt_metrics_to_yaml(df_metrics: pd.DataFrame, symprec: float) -> None:
    """Write geometry optimization metrics to model YAML metadata files.

    Args:
        df_metrics (pd.DataFrame): DataFrame with all geometry optimization metrics.
            Index = model names, columns = metric names including:
            - structure_rmsd_vs_dft: RMSD between predicted and DFT structures
            - n_sym_ops_mae: Mean absolute error in number of symmetry operations
            - symmetry_decrease: Fraction of structures with decreased symmetry
            - symmetry_match: Fraction of structures with matching symmetry
            - symmetry_increase: Fraction of structures with increased symmetry
            - n_structs: Number of structures evaluated
        symprec (float): symmetry precision for comparing ML and DFT relaxed structures.
    """
    for model_name in df_metrics.index:
        try:
            model = Model.from_label(model_name)
        except StopIteration:
            print(f"Skipping unknown {model_name=}")
            continue

        # Load existing metadata
        with open(model.yaml_path) as file:
            model_metadata = round_trip_yaml.load(file)

        all_metrics = model_metadata.setdefault("metrics", {})

        # Get metrics for this model
        model_metrics = df_metrics.loc[model_name]
        new_metrics = {
            str(Key.rmsd): float(round(model_metrics[MbdKey.structure_rmsd_vs_dft], 4)),
            str(Key.n_sym_ops_mae): float(round(model_metrics[Key.n_sym_ops_mae], 4)),
            str(Key.symmetry_decrease): float(
                round(model_metrics[Key.symmetry_decrease], 4)
            ),
            str(Key.symmetry_match): float(round(model_metrics[Key.symmetry_match], 4)),
            str(Key.symmetry_increase): float(
                round(model_metrics[Key.symmetry_increase], 4)
            ),
            str(Key.n_structs): int(model_metrics[Key.n_structs]),
        }
        symprec_key = f"{symprec=}".replace("e-0", "e-")

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
            Key.n_structs: "count",
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


def calc_geo_opt_metrics(df_geo_opt: pd.DataFrame) -> pd.DataFrame:
    """Calculate geometry optimization metrics for each model.

    Args:
        df_geo_opt (pd.DataFrame): DataFrame with geometry optimization metrics for all
            models and DFT reference. Must have a 2-level column MultiIndex with levels
            [model_name, property]. Required properties are:
            - structure_rmsd_vs_dft: RMSD between predicted and DFT structures
            - n_sym_ops_diff: Difference in number of symmetry operations vs DFT
            - spg_num_diff: Difference in space group number vs DFT

    Returns:
        pd.DataFrame: DataFrame with geometry optimization metrics. Shape = (n_models,
        n_metrics). Columns include:
        - structure_rmsd_vs_dft: Mean RMSD between predicted and DFT structures
        - n_sym_ops_mae: Mean absolute error in number of symmetry operations
        - symmetry_decrease: Fraction of structures with decreased symmetry
        - symmetry_match: Fraction of structures with matching symmetry
        - symmetry_increase: Fraction of structures with increased symmetry
        - n_structs: Number of structures evaluated
    """
    results: dict[str, dict[str, float]] = {}

    for model in set(df_geo_opt.columns.levels[0]) - {Key.dft.label}:
        try:
            # Get relevant columns for this model
            spg_diff = df_geo_opt[model][MbdKey.spg_num_diff]
            n_sym_ops_diff = df_geo_opt[model][MbdKey.n_sym_ops_diff]
            rmsd = df_geo_opt[model][MbdKey.structure_rmsd_vs_dft]

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

            results[model] = {
                str(MbdKey.structure_rmsd_vs_dft): float(mean_rmsd),
                str(Key.n_sym_ops_mae): float(sym_ops_mae),
                str(Key.symmetry_decrease): float(sym_decreased.sum() / total),
                str(Key.symmetry_match): float(sym_matched.sum() / total),
                str(Key.symmetry_increase): float(sym_increased.sum() / total),
                str(Key.n_structs): total,
            }
        except KeyError as exc:
            exc.add_note(
                f"Missing data for {model}, available columns={list(df_geo_opt[model])}"
            )
            raise

    return pd.DataFrame(results).T
