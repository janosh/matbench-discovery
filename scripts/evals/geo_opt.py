"""Generate geometry-optimization payloads from ML-vs-DFT structure analyses.

To compute analysis CSVs and write metrics for a new model, first run
scripts/analyze_geo_opt.py.

RMSD values (structure_rmsd_vs_dft) are normalized by (volume per atom)^1/3 and
thus unitless, not in Ångström. NaN values from structures that couldn't be
matched are filled with 1.0 (the stol value from StructureMatcher) for
metrics calculation.

TODO maybe add average_distance_within_threshold metric
https://github.com/FAIR-Chem/fairchem/blob/6329e922/src/fairchem/core/modules/evaluator.py#L318
"""

# %%
import numpy as np
import pandas as pd
import pymatviz as pmv
from pymatviz.enums import Key

from matbench_discovery import ROOT, figs, payload_numerics
from matbench_discovery.cli import cli_args, is_full_model_run
from matbench_discovery.data import (
    canonical_scientific_notation,
    file_ref_name,
    file_ref_url,
)
from matbench_discovery.enums import DataFiles, MbdKey, resolve_verified_file

symprec = 1e-5
metric_lvl = "metric"
# must match the key written by metrics/geo_opt.py, else lookups silently miss
symprec_str = f"symprec={canonical_scientific_notation(symprec)}"
df_dft_analysis = pd.read_csv(DataFiles.wbm_dft_geo_opt_symprec_1e_5.path, index_col=0)


# %% Load all model data
model_data: dict[str, pd.DataFrame] = {}
model_identities: dict[str, dict[str, object]] = {}
for model in cli_args.models:
    metrics = model.metrics.get("geo_opt") or {}
    if not file_ref_name(metrics.get("pred_file")):
        continue

    symprec_metrics = metrics.get(symprec_str, {})
    analysis_ref = symprec_metrics.get("analysis_file")
    if not (analysis_name := file_ref_name(analysis_ref)):
        print(f"Warning: {model.label} has no analysis file for {symprec_str}")
        continue

    analysis_path = f"{ROOT}/{analysis_name}"
    # fetch the (small) symmetry-analysis CSV from figshare if missing so payloads can
    # be regenerated on machines/CI that never ran the model locally
    resolve_verified_file(
        analysis_path,
        url=file_ref_url(analysis_ref),
        label=model.label,
        md5=analysis_ref.get("md5") if isinstance(analysis_ref, dict) else None,
    )

    rel_path = analysis_path.rsplit("/models/", maxsplit=1)[-1]
    print(f"Found {model.label} analysis file:\n  {rel_path}")
    model_data[model.key] = pd.read_csv(analysis_path, index_col=0)
    model_identities[model.key] = figs.model_payload_identity(
        model, "geo_opt_analysis", analysis_path
    )

if not model_data:
    print("No geo-opt analysis data for requested models, nothing to do")
    # on subset runs (e.g. ingesting an energy-only model) there's nothing to merge
    # into the site payloads, which is fine; a full run with no data is a config error
    raise SystemExit(1 if is_full_model_run() else 0)

model_keys = tuple(model_data)
print(f"\nLoaded {len(model_data)=} models, joined with DFT data into df_all")
df_all = pd.concat(
    model_data | {Key.dft.label: df_dft_analysis},
    axis="columns",
    names=["model", "metric"],
)
del model_data


# %% ML-vs-DFT relaxed spacegroup correspondence
df_spg = df_all.xs(Key.spg_num, level=metric_lvl, axis="columns").convert_dtypes()
spg_sankey_models: list[dict[str, object]] = []
for model_key in model_keys:
    # get most common pairs of DFT/Model spacegroups
    common_dft_spgs, common_model_spgs = zip(
        *df_spg[[Key.dft.label, model_key]].value_counts().head(10).index,
        strict=True,
    )
    df_spg_common = df_spg[
        df_spg[Key.dft.label].isin(common_dft_spgs)
        & df_spg[model_key].isin(common_model_spgs)
    ].sort_values(by=Key.dft.label)
    flow_data = pmv.sankey_flow_data(
        df_spg_common.reset_index(),
        [Key.dft.label, model_key],
    )
    spg_sankey_models.append(
        {
            **model_identities[model_key],
            **payload_numerics.sankey_payload_from_flow(flow_data),
        }
    )


def payload_provenance(
    parameters: dict[str, object], *, uses_payload_numerics: bool = False
) -> dict[str, object]:
    """Build shared provenance for one geo-opt payload."""
    return figs.build_payload_provenance(
        generator=__file__,
        benchmark_inputs={
            "wbm_dft_geo_opt_symprec_1e_5": (
                DataFiles.wbm_dft_geo_opt_symprec_1e_5.path
            )
        },
        source_files=(
            {"payload_numerics": payload_numerics.__file__}
            if uses_payload_numerics
            else {}
        ),
        parameters={"symprec": symprec} | parameters,
        packages=("numpy", "pandas", "pymatviz"),
    )


# one combined payload for all models (writer sorts entries for determinism)
figs.write_site_payload(
    "spg-sankeys",
    {
        **payload_provenance(
            {
                "coordinate_decimals": payload_numerics.COORD_DECIMALS,
                "top_spacegroup_pairs": 10,
            },
            uses_payload_numerics=True,
        ),
        "models": spg_sankey_models,
    },
)


# %% Difference in symmetry-operation counts for ML vs DFT-relaxed structures
df_sym_ops_diff = df_all.drop(Key.dft.label, level=Key.model, axis="columns").xs(
    MbdKey.n_sym_ops_diff, level=metric_lvl, axis="columns"
)
models_by_std = sorted(df_sym_ops_diff.std().items(), key=lambda item: item[1])
sym_ops_models: list[dict[str, object]] = []
for raw_model_key, std_dev in models_by_std:
    # sort_index + int cast: multi-model runs concat to a NaN-padded float column but
    # single-model runs stay int, which flips both value_counts' tie order and the
    # 0.0-vs-0 JSON repr, breaking merge == full-regen byte identity
    model_key = str(raw_model_key)
    value_counts = df_sym_ops_diff[model_key].value_counts().sort_index()
    sym_ops_models.append(
        {
            **model_identities[model_key],
            "sigma": float(f"{std_dev:.3g}"),
            "x": [int(val) for val in value_counts.index],
            "y": value_counts.tolist(),
        }
    )

figs.write_site_payload(
    "sym-ops-diff-bar",
    {
        **payload_provenance({"sigma_significant_digits": 3}),
        "models": sym_ops_models,
    },
)


# %% Cumulative distribution of RMSD vs DFT
x_max = 0.05
rmsd_cdf_models: list[dict[str, object]] = []
for model_key in model_keys:
    rmsd_vals = df_all.xs(
        MbdKey.structure_rmsd_vs_dft, level=metric_lvl, axis="columns"
    )[model_key].dropna()

    # Calculate CDF
    sorted_rmsds = np.sort(rmsd_vals)
    cumulative = np.arange(1, len(sorted_rmsds) + 1) / len(sorted_rmsds)

    # Sample at most 200 distinct points evenly across the range
    indices = payload_numerics.evenly_spaced_indices(len(sorted_rmsds), 200)
    sampled_rmsds = sorted_rmsds[indices]
    sampled_cumulative = cumulative[indices]

    # calculate AUC only up to x_max
    auc = (
        np.trapezoid(
            sampled_cumulative[sampled_rmsds <= x_max],
            sampled_rmsds[sampled_rmsds <= x_max],
        )
        / x_max
    )
    rmsd_cdf_models.append(
        {
            **model_identities[model_key],
            "auc": float(f"{auc:.3g}"),
            "x": payload_numerics.round_list(sampled_rmsds),
            "y": payload_numerics.round_list(sampled_cumulative),
        }
    )

figs.write_site_payload(
    "struct-rmsd-cdf",
    {
        **payload_provenance(
            {
                "auc_significant_digits": 3,
                "cdf_samples": 200,
                "coordinate_decimals": payload_numerics.COORD_DECIMALS,
                "x_max": x_max,
            },
            uses_payload_numerics=True,
        ),
        "models": rmsd_cdf_models,
    },
)
