"""Functions to calculate and save geometry optimization metrics."""

import pandas as pd
from pymatviz.enums import Key
from ruamel.yaml.comments import CommentedMap

from matbench_discovery import ROOT
from matbench_discovery.data import update_yaml_at_path
from matbench_discovery.enums import MbdKey, Model


def write_metrics_to_yaml(
    metrics: dict[str, float] | pd.DataFrame,
    model: Model,
    analysis_file_path: str,
    section: str = "symmetry",
    symprec: float | None = None,
) -> dict[str, str | float]:
    """Write metrics to model's YAML file under the specified section.

    Args:
        metrics (dict[str, float] | pd.DataFrame): Either a DataFrame with geometry
            optimization metrics columns or a dictionary of metrics
        model (Model): Model to write metrics for
        analysis_file_path (str): Path to analysis file (should be the combined CSV)
        section (str): Section name to write metrics to. Either "symmetry" or "distance"
        symprec (float, optional): Symmetry precision used for analysis.
            Required when section="symmetry", ignored when section="distance".

    Returns:
        dict[str, str | float]: Metrics added to the YAML for this model
    """
    # Convert absolute path to relative path if needed
    analysis_file_path = analysis_file_path.removeprefix(f"{ROOT}/")

    # Convert DataFrame to dict if needed
    if isinstance(metrics, pd.DataFrame):
        metrics = calc_geo_opt_metrics(metrics)

    # Only include relevant metrics for the specified section
    formatted_metrics = CommentedMap()

    # Define units for different metrics
    metric_units = {
        str(Key.rmsd): "unitless",
        str(Key.n_sym_ops_mae): "unitless",
        str(Key.symmetry_decrease): "fraction",
        str(Key.symmetry_match): "fraction",
        str(Key.symmetry_increase): "fraction",
        str(Key.n_structures): "count",
    }

    # Ensure common metadata exists at the root level
    root = update_yaml_at_path(model.yaml_path, "metrics.geo_opt")
    if "analysis_file" not in root:
        root["analysis_file"] = analysis_file_path
        root["analysis_file_url"] = None  # to be filled after upload to figshare
        update_yaml_at_path(model.yaml_path, "metrics.geo_opt", root)

    # Handle n_structures at top level if present in metrics
    if str(Key.n_structures) in metrics and str(Key.n_structures) not in root:
        root[str(Key.n_structures)] = int(metrics[str(Key.n_structures)])
        update_yaml_at_path(model.yaml_path, "metrics.geo_opt", root)
        # Add unit comment if possible
        yaml_root = update_yaml_at_path(model.yaml_path, "metrics.geo_opt")
        if hasattr(yaml_root, "yaml_add_eol_comment"):
            yaml_root.yaml_add_eol_comment("count", str(Key.n_structures), column=1)

    # Build section-specific metrics
    for key, value in metrics.items():
        if (
            section == "symmetry"
            and key != str(MbdKey.structure_rmsd_vs_dft)
            and key != str(Key.n_structures)
        ):
            if isinstance(value, float):
                formatted_metrics[key] = round(value, 4)
            else:
                formatted_metrics[key] = value
        elif section == "distance" and key == str(MbdKey.structure_rmsd_vs_dft):
            formatted_metrics[str(Key.rmsd)] = round(float(value), 4)

    # Add units as YAML end-of-line comments
    for key in formatted_metrics:
        if key in metric_units and hasattr(formatted_metrics, "yaml_add_eol_comment"):
            formatted_metrics.yaml_add_eol_comment(metric_units[key], key, column=1)

    # Update the appropriate section in the YAML
    if section == "symmetry":
        if symprec is None:
            raise ValueError("symprec must be provided when section='symmetry'")
        symprec_key = f"{symprec=:.0e}".replace("e-0", "e-")

        # Ensure symmetry section exists
        update_yaml_at_path(model.yaml_path, "metrics.geo_opt.symmetry", {})

        update_yaml_at_path(
            model.yaml_path,
            f"metrics.geo_opt.{section}.{symprec_key}",
            formatted_metrics,
        )
    else:
        # Ensure distance section exists
        update_yaml_at_path(model.yaml_path, "metrics.geo_opt.distance", {})

        update_yaml_at_path(
            model.yaml_path, f"metrics.geo_opt.{section}", formatted_metrics
        )

    return formatted_metrics


def calc_geo_opt_metrics(df_model_analysis: pd.DataFrame) -> dict[str, float]:
    """Calculate geometry optimization metrics for a single model.

    Args:
        df_model_analysis (pd.DataFrame): DataFrame with geometry optimization metrics
            for one model. Depending on the analysis type, it may contain columns:
            - structure_rmsd_vs_dft: RMSD between predicted and DFT structures
            - n_sym_ops_diff: Difference in number of symmetry operations vs DFT
            - spg_num_diff: Difference in space group number vs DFT

    Returns:
        dict[str, float]: Geometry optimization metrics including:
            - structure_rmsd_vs_dft: Mean RMSD between predicted and DFT structures
            - n_sym_ops_mae: Mean absolute error in number of symmetry operations
            - symmetry_decrease/match/increase: Fractions of total structures
            - n_structures: Number of structures evaluated
    """
    result: dict[str, float] = {str(Key.n_structures): float(len(df_model_analysis))}

    # Calculate RMSD metrics if distance analysis was performed
    if MbdKey.structure_rmsd_vs_dft in df_model_analysis:
        # First convert to numeric type, then fill NaN values
        rmsd_vals = pd.to_numeric(
            df_model_analysis[MbdKey.structure_rmsd_vs_dft], errors="coerce"
        )
        # Fill NaN values with 1.0 (the stol value we set in StructureMatcher)
        mean_rmsd = rmsd_vals.fillna(1.0).mean()
        result[str(MbdKey.structure_rmsd_vs_dft)] = float(mean_rmsd)

    # Calculate symmetry metrics if columns exist
    has_symmetry = (
        MbdKey.spg_num_diff in df_model_analysis
        and MbdKey.n_sym_ops_diff in df_model_analysis
    )

    if has_symmetry:
        spg_diff = df_model_analysis[MbdKey.spg_num_diff]
        n_sym_ops_diff = df_model_analysis[MbdKey.n_sym_ops_diff]

        # Get valid count for percentage calculations
        valid_count = float(len(spg_diff.dropna()))

        if valid_count > 0:
            # Calculate MAE and symmetry metrics
            result[str(Key.n_sym_ops_mae)] = float(n_sym_ops_diff.abs().mean())

            # Calculate symmetry change categories
            changed_mask = spg_diff != 0
            sym_decreased = (n_sym_ops_diff < 0) & changed_mask
            sym_increased = (n_sym_ops_diff > 0) & changed_mask

            # Calculate fractions
            result[str(Key.symmetry_decrease)] = float(
                sym_decreased.sum() / valid_count
            )
            result[str(Key.symmetry_match)] = float((~changed_mask).sum() / valid_count)
            result[str(Key.symmetry_increase)] = float(
                sym_increased.sum() / valid_count
            )

    return result
