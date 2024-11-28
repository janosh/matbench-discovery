"""Functions to calculate and save geometry optimization metrics."""

import pandas as pd
from pymatviz.enums import Key, Task
from ruamel.yaml.comments import CommentedMap

from matbench_discovery.data import Model, round_trip_yaml
from matbench_discovery.enums import MbdKey


def write_geo_opt_metrics_to_yaml(df_metrics: pd.DataFrame) -> None:
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

        geo_opt_metrics = CommentedMap(
            all_metrics.setdefault(Task.geo_opt, {}) | new_metrics
        )

        # Define units for metrics
        metric_units = {
            Key.rmsd: "Å",
            Key.n_sym_ops_mae: "count",
            Key.symmetry_decrease: "fraction",
            Key.symmetry_match: "fraction",
            Key.symmetry_increase: "fraction",
            Key.n_structs: "count",
        }

        # Add units as YAML end-of-line comments
        for key in new_metrics:
            if unit := metric_units.get(key):
                geo_opt_metrics.yaml_add_eol_comment(unit, key, column=1)

        all_metrics[Task.geo_opt] = geo_opt_metrics

        # Write back to file
        with open(model.yaml_path, mode="w") as file:
            round_trip_yaml.dump(model_metadata, file)


def analyze_symmetry_changes(df_sym: pd.DataFrame) -> pd.DataFrame:
    """Analyze how often each model's predicted structure has different symmetry vs DFT.

    Returns:
        pd.DataFrame: DataFrame with columns for fraction of structures where symmetry
            decreased, matched, or increased vs DFT.
    """
    results: dict[str, dict[str, float]] = {}

    for model in df_sym.columns.levels[0]:
        if model == Key.dft.label:  # don't compare DFT to itself
            continue
        try:
            spg_diff = df_sym[model][MbdKey.spg_num_diff]
            n_sym_ops_diff = df_sym[model][MbdKey.n_sym_ops_diff]
            total = len(spg_diff.dropna())

            # Count cases where spacegroup changed
            changed_mask = spg_diff != 0
            # Among changed cases, count whether symmetry increased or decreased
            sym_decreased = (n_sym_ops_diff < 0) & changed_mask
            sym_increased = (n_sym_ops_diff > 0) & changed_mask
            sym_matched = ~changed_mask

            results[model] = {
                str(Key.symmetry_decrease): float(sym_decreased.sum() / total),
                str(Key.symmetry_match): float(sym_matched.sum() / total),
                str(Key.symmetry_increase): float(sym_increased.sum() / total),
                str(Key.n_structs): total,
            }
        except KeyError as exc:
            exc.add_note(
                f"Missing data for {model}, available columns={list(df_sym[model])}"
            )
            raise

    return pd.DataFrame(results).T
