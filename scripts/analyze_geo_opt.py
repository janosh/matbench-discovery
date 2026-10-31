"""Run this script to add/update geometry optimization analysis for new models to
individual CSV files (one per symprec value) in a model's directory.

Output files will have the same name as the input file containing the model's relaxed
structures, but with the symprec value appended to the filename.

Note: StructureMatcher is configured with stol=1.0 and scale=False for exact matching.
RMSD values are normalized (unitless) and NaN values are filled with 1.0 (the stol
value) to properly account for structures that couldn't be matched.

Example usage (moyopy comes from the symmetry extra):
    uv run --with-editable '.[symmetry]' scripts/analyze_geo_opt.py \
        --models mace_mp_0 m3gnet --symprec 1e-2 1e-5
    # only analyze first 10 structures (CSV goes to a debug/ subdir, no YAML write)
    uv run --with-editable '.[symmetry]' scripts/analyze_geo_opt.py --debug 10
    uv run --with-editable '.[symmetry]' scripts/analyze_geo_opt.py --workers 4
"""

# %%
import importlib
import importlib.metadata
import itertools
import os
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from typing import TYPE_CHECKING, Any, Final

import pandas as pd
from pymatgen.core import Structure
from pymatviz.enums import Key

from matbench_discovery import ROOT, file_digest
from matbench_discovery.cli import cli_parser, parse_cli_args
from matbench_discovery.data import (
    artifact_filename,
    canonical_scientific_notation,
    file_ref_name,
    file_ref_url,
    parse_artifact_filename,
)
from matbench_discovery.enums import DataFiles, Model
from matbench_discovery.metrics import geo_opt
from matbench_discovery.remote.fetch import maybe_auto_download_file
from matbench_discovery.structure import symmetry

if TYPE_CHECKING:
    from collections.abc import Sequence

DFT_REFERENCE_FILES: Final[dict[float, DataFiles]] = {
    1e-2: DataFiles.wbm_dft_geo_opt_symprec_1e_2,
    1e-5: DataFiles.wbm_dft_geo_opt_symprec_1e_5,
}


def load_model_structs(model: Model, debug_mode: int) -> dict[str, Structure]:
    """Download (if needed), verify and parse a model's geo-opt JSON Lines into
    Structures keyed by ID (slow, GBs of RAM): only for uncached analyses.
    """
    geo_opt_path = model.geo_opt_path
    if geo_opt_path is None:
        raise ValueError(f"{model.label} has no metrics.geo_opt.pred_file")
    try:
        df_ml_structs = pd.read_json(geo_opt_path, lines=True)
    except Exception as exc:
        exc.add_note(f"{model.label=} {geo_opt_path=}")
        raise

    if missing_columns := {"material_id", "structure"} - set(df_ml_structs):
        raise ValueError(
            f"Missing canonical geo-opt columns {sorted(missing_columns)} in "
            f"{geo_opt_path}"
        )
    material_ids = df_ml_structs["material_id"]
    if material_ids.isna().any():
        raise ValueError(f"Null material_id values in {geo_opt_path}")
    if dupes := sorted(material_ids[material_ids.duplicated()].astype(str).unique()):
        raise ValueError(f"Duplicate material_id values in {geo_opt_path}: {dupes}")
    df_ml_structs = df_ml_structs.set_index("material_id")

    if debug_mode:
        df_ml_structs = df_ml_structs.head(debug_mode)

    return {
        mat_id: Structure.from_dict(struct_dict)
        for mat_id, struct_dict in df_ml_structs["structure"].items()
    }


def analyze_model_symprec(
    model: Model,
    symprec: float,
    moyo_version: str,
    df_dft_analysis: pd.DataFrame,
    dft_structs: dict[str, Structure],
    *,
    debug_mode: int = 0,
    pbar_pos: int = 0,  # tqdm progress bar position
    overwrite: bool = False,
) -> pd.DataFrame | None:
    """Analyze a single model for a single symprec value."""
    geo_opt_metrics: dict[str, Any] = model.metrics.get("geo_opt") or {}

    # only the registered name: a cached analysis never needs the multi-GB file
    geo_opt_name = file_ref_name(geo_opt_metrics.get("pred_file"))
    if not geo_opt_name:
        print(f"⚠️ {model.label} has no relaxed structures file")
        return None
    geo_opt_path = f"{ROOT}/{geo_opt_name}"

    symprec_str = f"symprec={canonical_scientific_notation(symprec)}"
    geo_opt_basename = os.path.basename(geo_opt_path)
    if parse_artifact_filename(geo_opt_basename) != "geo_opt":
        raise ValueError(f"Expected geo-opt artifact, got {geo_opt_path!r}")
    analysis_name = artifact_filename(
        geo_opt_basename[:10],
        "geo_opt_analysis",
        symprec=symprec,
        moyo_version=moyo_version,
    )
    geo_opt_dir = os.path.dirname(geo_opt_path)
    if debug_mode:
        # a debug subset must never clobber the full analysis CSV (nor, below, the
        # model YAML), so it lands in a sibling debug/ dir like run_discovery's dry-run/
        # keyed by subset size so --debug 3 doesn't reuse a cached --debug 10 CSV
        geo_opt_csv_path = f"{geo_opt_dir}/debug/{debug_mode}/{analysis_name}"
    else:
        geo_opt_csv_path = f"{geo_opt_dir}/{analysis_name}"

        # only reuse/download a published analysis file if its path matches exactly
        # `or {}`: a null YAML value (symprec=1e-2: ~) must read as no metrics
        symprec_metrics = geo_opt_metrics.get(symprec_str) or {}
        analysis_ref = symprec_metrics.get("analysis_file")
        analysis_url = file_ref_url(analysis_ref)
        analysis_file = file_ref_name(analysis_ref)
        analysis_file_path = f"{ROOT}/{analysis_file}" if analysis_file else ""

        if analysis_file_path == geo_opt_csv_path:
            if not os.path.isfile(geo_opt_csv_path) and analysis_url:
                maybe_auto_download_file(
                    analysis_url,
                    geo_opt_csv_path,
                    label=f"{model.label} {symprec_str}",
                )
        elif analysis_file:
            # Paths differ (moyo version change, etc.) - will recompute
            print(f"⚠️ {model.label} {symprec_str=} path mismatch, will recompute")
            print(f"  - Expected: {geo_opt_csv_path}")
            print(f"  - Found: {analysis_file_path}")

    if os.path.isfile(geo_opt_csv_path) and not overwrite:
        # reuse the cached analysis but still (re)derive metrics from it below so the
        # YAML never lags behind an already-analyzed CSV
        print(f"{model.label} already analyzed at {geo_opt_csv_path}")
        df_ml_geo_analysis = pd.read_csv(geo_opt_csv_path, index_col=0)
    else:
        # only reached without a cached CSV or with overwrite=True, so an existing
        # file here is always being overwritten
        action = "Overwriting" if os.path.isfile(geo_opt_csv_path) else "Analyzing"
        print(f"{action} {model.label} for {symprec=}")
        model_structs = load_model_structs(model, debug_mode)

        pbar_desc = f"Process {pbar_pos}: Analyzing {model.label} for {symprec=}"
        df_model_analysis = symmetry.get_sym_info_from_structs(
            model_structs,
            pbar=dict(desc=pbar_desc, position=pbar_pos, leave=True),
            symprec=symprec,
        )

        pbar_desc = f"Process {pbar_pos}: Comparing DFT vs {model.label} for {symprec=}"
        df_ml_geo_analysis = symmetry.pred_vs_ref_struct_symmetry(
            df_model_analysis,
            df_dft_analysis,
            model_structs,
            dft_structs,
            pbar=dict(desc=pbar_desc, position=pbar_pos, leave=True),
        )

        os.makedirs(os.path.dirname(geo_opt_csv_path), exist_ok=True)
        df_ml_geo_analysis.to_csv(geo_opt_csv_path)
        print(f"Completed {model.label} {symprec=}, saved to {geo_opt_csv_path}")

    metrics_dict = geo_opt.calc_geo_opt_metrics(df_ml_geo_analysis)
    df_metrics = pd.DataFrame([metrics_dict])
    if debug_mode:
        print(
            f"Debug mode ({debug_mode} structures): not writing metrics.geo_opt for "
            f"{model.label} {symprec=} to YAML\n{df_metrics.T}"
        )
    else:
        geo_opt.write_metrics_to_yaml(df_metrics, model, symprec, geo_opt_csv_path)
    return df_metrics


def write_dft_reference(
    dft_structs: dict[str, Structure], symprec: float, moyo_version: str, out_dir: str
) -> pd.DataFrame:
    """Recompute the WBM DFT symmetry reference into a new file named by date and
    moyopy version (never overwriting one) to upload and register in data-files.yml.
    """
    out_path = (
        f"{out_dir}/wbm/{datetime.now(tz=UTC).date()}-dft-geo-opt-symprec="
        f"{canonical_scientific_notation(symprec)}-moyo={moyo_version}.csv.gz"
    )
    if os.path.isfile(out_path):
        raise FileExistsError(f"Refusing to overwrite DFT reference {out_path}")
    df_dft_analysis = symmetry.get_sym_info_from_structs(
        dft_structs,
        pbar=dict(desc=f"Getting DFT symmetries {symprec=}"),
        symprec=symprec,
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_dft_analysis.to_csv(out_path)
    md5 = file_digest(out_path, "md5")
    print(f"Wrote {out_path} ({md5=}): upload and register it in data-files.yml")
    return df_dft_analysis


def main(raw_args: list[str] | None = None) -> int:
    """Analyze geometry optimization results.

    Args:
        raw_args: Command line arguments. Defaults to sys.argv[1:].

    Returns:
        Exit code (0 for success).
    """
    geo_opt_group = cli_parser.add_argument_group(
        "geo_opt", "Arguments for geometry optimization analysis"
    )
    geo_opt_group.add_argument(
        "--symprec",
        nargs="+",
        type=float,
        default=[1e-2, 1e-5],
        help="Symmetry precision values to analyze.",
    )
    geo_opt_group.add_argument(
        "--write-dft-reference",
        action="store_true",
        help="Recompute the WBM DFT symmetry reference with the installed moyopy.",
    )
    args = parse_cli_args(raw_args)

    # > 0 activates debug mode, analyzing only that many structures
    debug_mode: Final[int] = args.debug
    symprec_values: Final[Sequence[float]] = args.symprec
    moyo_version = importlib.metadata.version("moyopy")

    # %%
    print("Loading WBM PBE structures...")
    wbm_cse_path = DataFiles.wbm_computed_structure_entries.path
    df_wbm_structs: pd.DataFrame = pd.read_json(
        wbm_cse_path, lines=True, orient="records"
    ).set_index(Key.mat_id)

    dft_structs: dict[str, Structure] = {
        mat_id: Structure.from_dict(cse[Key.structure])
        for mat_id, cse in df_wbm_structs[Key.computed_structure_entry].items()
    }

    # %% Load the registered DFT reference for each symprec value. Model analyses are
    # only comparable to a reference computed with the same moyopy version, so a version
    # mismatch fails here instead of silently recomputing an unregistered reference.
    dft_analysis_dict: dict[float, pd.DataFrame] = {}
    for symprec in symprec_values:
        symprec_str = f"symprec={canonical_scientific_notation(symprec)}"
        data_file = DFT_REFERENCE_FILES.get(symprec)
        if data_file is None:
            raise ValueError(
                f"No registered DFT reference for {symprec_str}, available: "
                f"{sorted(DFT_REFERENCE_FILES)}"
            )
        if args.write_dft_reference:
            dft_analysis_dict[symprec] = write_dft_reference(
                dft_structs, symprec, moyo_version, DataFiles.base_dir
            )
            continue
        if f"moyo={moyo_version}" not in data_file.rel_path:
            raise ValueError(
                f"Installed moyopy {moyo_version} does not match the registered "
                f"DFT reference {data_file.rel_path}. Install the registered version, "
                "or pass --write-dft-reference and register the new file in "
                "data-files.yml"
            )
        dft_analysis_dict[symprec] = pd.read_csv(data_file.path).set_index(Key.mat_id)

    tasks = list(itertools.product(args.models, symprec_values))
    n_workers = min(len(tasks), args.workers)

    # %% Process model-symprec combinations in parallel
    print(
        f"\nAnalyzing {len(tasks)} model-symprec combos using {n_workers} processes..."
    )

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(
                analyze_model_symprec,
                model=model_name,
                symprec=symprec,
                moyo_version=moyo_version,
                df_dft_analysis=dft_analysis_dict[symprec],
                dft_structs=dft_structs,
                debug_mode=debug_mode,
                pbar_pos=idx,  # assign unique position to each task's progress bar
                overwrite=args.overwrite,
            )
            for idx, (model_name, symprec) in enumerate(tasks)
        ]
        for future in futures:
            future.result()  # re-raises any worker exception

    print(f"\nAll {len(tasks)} model-symprec combinations processed!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
