"""Functions to calculate and save geometry optimization metrics."""

import numpy as np
import pandas as pd
from pymatviz.enums import Key, Task
from ruamel.yaml.comments import CommentedMap

from matbench_discovery.data import Model, round_trip_yaml
from matbench_discovery.enums import MbdKey


def write_geo_opt_metrics_to_yaml(
    df_sym: pd.DataFrame, df_sym_changes: pd.DataFrame
) -> None:
    """Write geometry optimization metrics to model YAML metadata files."""
    for model_name in df_sym.columns.levels[0]:
        try:
            model = Model.from_label(model_name)
        except StopIteration:
            print(f"Skipping {model_name}")
            continue

        df_rmsd = df_sym.xs(
            MbdKey.structure_rmsd_vs_dft, level=MbdKey.sym_prop, axis="columns"
        ).round(4)
        if model.label not in df_rmsd:
            print(f"No RMSD column for {model.label}")
            return

        # Calculate RMSD
        rmsd = df_rmsd[model.label].mean(axis=0)

        # Add symmetric mean error calculation for number of symmetry operations
        sym_ops_diff = df_sym.drop(Key.dft.label, level=Key.model, axis="columns")[
            model.label
        ][MbdKey.n_sym_ops_diff]

        # Calculate symmetric mean error for each model
        sym_ops_mae = np.mean(np.abs(sym_ops_diff))

        # Calculate symmetry change statistics
        if model.label not in df_sym_changes.index:
            print(f"No symmetry data for {model.label}")
            return
        sym_changes = df_sym_changes.round(4).loc[model.label].to_dict()
        # Combine metrics
        with open(model.yaml_path) as file:  # Load existing metadata
            model_metadata = round_trip_yaml.load(file)

        all_metrics = model_metadata.setdefault("metrics", {})
        new_metrics = {
            str(Key.rmsd): float(round(rmsd, 4)),
            str(Key.n_sym_ops_mae): float(round(sym_ops_mae, 4)),
            **sym_changes,
        }
        geo_opt_metrics = CommentedMap(
            all_metrics.setdefault(Task.geo_opt, {}) | new_metrics
        )
        metric_units = dict.fromkeys(sym_changes, "fraction") | {
            Key.rmsd: "Å",
            Key.n_structs: "count",
        }
        # Add units as YAML end-of-line comments
        for key in new_metrics:
            if unit := metric_units.get(key):
                geo_opt_metrics.yaml_add_eol_comment(unit, key, column=1)

        all_metrics[Task.geo_opt] = geo_opt_metrics

        with open(model.yaml_path, mode="w") as file:  # Write back to file
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
