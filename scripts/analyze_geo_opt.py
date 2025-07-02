"""Run this script to add/update geometry optimization analysis for new models to
individual CSV files (one per symprec value) in a model's directory.

Output files will have the same name as the input file containing the model's relaxed
structures, but with the symprec value appended to the filename.

Note: StructureMatcher is configured with stol=1.0 and scale=False for exact matching.
RMSD values are normalized (unitless) and NaN values are filled with 1.0 (the stol
value) to properly account for structures that couldn't be matched.

Example usage:
    python scripts/analyze_geo_opt.py --models mace_mp_0 m3gnet --symprec 1e-2 1e-5
    python scripts/analyze_geo_opt.py --debug 10  # only analyze first 10 structures
    python scripts/analyze_geo_opt.py --workers 4  # use 4 CPU cores
"""

# %%
import importlib
import importlib.metadata
import itertools
import os
from concurrent.futures import ProcessPoolExecutor
from typing import TYPE_CHECKING, Any, Final

import pandas as pd
from pymatgen.core import Structure
from pymatviz.enums import Key

from matbench_discovery import ROOT
from matbench_discovery.cli import cli_parser
from matbench_discovery.data import update_yaml_file
from matbench_discovery.enums import DataFiles, Model
from matbench_discovery.metrics import geo_opt
from matbench_discovery.structure import symmetry

if TYPE_CHECKING:
    from collections.abc import Sequence


def analyze_model_symprec(
    model: Model,
    symprec: float,
    moyo_version: str,
    df_dft_analysis: pd.DataFrame,
    dft_structs: dict[str, Structure],
    *,
    debug_mode: int = 0,
    pbar_pos: int = 0,  # tqdm progress bar position
    overwrite: bool = False,  # Whether to overwrite existing analysis files
) -> pd.DataFrame | None:
    """Analyze a single model for a single symprec value."""
    geo_opt_metrics: dict[str, Any] = model.metadata.get("metrics", {}).get(
        "geo_opt", {}
    )

    # skip models that don't support geometry optimization
    if geo_opt_metrics in ("not applicable", "not available"):
        print(f"⚠️ {model.label} does not support geometry optimization")
        return None

    if not model.geo_opt_path:
        print(f"⚠️ {model.label} has no relaxed structures file")
        return None

    if not os.path.isfile(ml_relaxed_structs_path := model.geo_opt_path):
        print(
            f"⚠️ {model.label}-relaxed structures not found, expected "
            f"at {ml_relaxed_structs_path}"
        )
        return None

    # Convert to JSON lines format if needed
    jsonl_path = ml_relaxed_structs_path
    if not ml_relaxed_structs_path.endswith((".jsonl", ".jsonl.gz")):
        # Try reading as JSON lines first (file might already be in correct format)
        try:
            pd.read_json(ml_relaxed_structs_path, lines=True)
        except ValueError:
            # Not JSON lines format, convert it
            jsonl_path = ml_relaxed_structs_path.rsplit(".", 2)[0] + ".jsonl.gz"
            if not os.path.isfile(jsonl_path):
                print(f"Converting {ml_relaxed_structs_path} to JSON lines format...")
                df = pd.read_json(ml_relaxed_structs_path)
                df.to_json(jsonl_path, orient="records", lines=True, compression="gzip")
                # Update model yaml with new path
                yaml_path = f"models/{model.name}/{model.name}.yml"
                if os.path.isfile(yaml_path):
                    update_yaml_file(
                        yaml_path,
                        "metrics.geo_opt.struct_file",
                        {"metrics": {"geo_opt": {"struct_file": jsonl_path}}},
                    )
            ml_relaxed_structs_path = jsonl_path

    # Load model structures
    try:
        df_ml_structs = pd.read_json(ml_relaxed_structs_path, lines=True)
    except Exception as exc:
        exc.add_note(f"{model.label=} {ml_relaxed_structs_path=}")
        raise

    # try normalize material ID column or raise
    if Key.mat_id in df_ml_structs:
        df_ml_structs = df_ml_structs.set_index(Key.mat_id)
    elif df_ml_structs.index[0].startswith("wbm-"):
        df_ml_structs.index.name = Key.mat_id
        df_ml_structs.reset_index().to_json(ml_relaxed_structs_path)
    else:
        raise ValueError(f"Could not infer ID column from {df_ml_structs.columns}")

    if debug_mode:
        df_ml_structs = df_ml_structs.head(debug_mode)

    struct_col = geo_opt_metrics.get("struct_col")
    if struct_col not in df_ml_structs:
        struct_cols = [col for col in df_ml_structs if Key.structure in col]
        print(
            f"⚠️ {struct_col=} not found in {model.label}-relaxed structures loaded "
            f"from {ml_relaxed_structs_path}. Did you mean one of {struct_cols}?"
        )
        return None

    # Convert structures
    model_structs = {
        mat_id: Structure.from_dict(struct_dict)
        for mat_id, struct_dict in df_ml_structs[struct_col].items()
    }

    symprec_str = f"symprec={symprec:.0e}".replace("e-0", "e-")
    # Remove common file extensions properly
    geo_opt_filename = model.geo_opt_path
    for suffix in [".jsonl.gz", ".json.gz", ".jsonl", ".json"]:
        if geo_opt_filename.endswith(suffix):
            geo_opt_filename = geo_opt_filename.removesuffix(suffix)
            break
    geo_opt_csv_path = f"{geo_opt_filename}-{symprec_str}-{moyo_version}.csv.gz"

    if os.path.isfile(geo_opt_csv_path) and not overwrite:
        print(f"{model.label} already analyzed at {geo_opt_csv_path}")
        return pd.read_csv(geo_opt_csv_path)

    action = (
        "Overwriting" if overwrite and os.path.isfile(geo_opt_csv_path) else "Analyzing"
    )
    print(f"{action} {model.label} for {symprec=}")

    # Analyze symmetry for current symprec
    pbar_desc = f"Process {pbar_pos}: Analyzing {model.label} for {symprec=}"
    df_model_analysis = symmetry.get_sym_info_from_structs(
        model_structs,
        pbar=dict(desc=pbar_desc, position=pbar_pos, leave=True),
        symprec=symprec,
    )

    # Compare with DFT reference
    pbar_desc = f"Process {pbar_pos}:Comparing DFT vs {model.label} for {symprec=}"
    # break here
    df_ml_geo_analysis = symmetry.pred_vs_ref_struct_symmetry(
        df_model_analysis,
        df_dft_analysis,
        model_structs,
        dft_structs,
        pbar=dict(desc=pbar_desc, position=pbar_pos, leave=True),
    )

    # Save model results
    df_ml_geo_analysis.to_csv(geo_opt_csv_path)
    print(f"Completed {model.label} {symprec=} and saved results to {geo_opt_csv_path}")

    # Calculate metrics and write to YAML
    df_metrics = geo_opt.calc_geo_opt_metrics(df_ml_geo_analysis)
    geo_opt.write_metrics_to_yaml(df_metrics, model, symprec, geo_opt_csv_path)
    return df_metrics


def main(args: list[str] | None = None) -> int:
    """Main function to analyze geometry optimization results.

    Args:
        args: Command line arguments. If None, sys.argv[1:] will be used.

    Returns:
        int: Exit code (0 for success).
    """
    # Add geo_opt specific arguments to the central CLI parser
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
    args, _unknown = cli_parser.parse_known_args(args)

    # set to > 0 to activate debug mode, only that many structures will be analyzed
    debug_mode: Final[int] = args.debug
    # List of symprec values to analyze
    symprec_values: Final[Sequence[float]] = args.symprec

    # Get list of models to analyze
    moyo_version = f"moyo={importlib.metadata.version('moyopy')}"

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

    # %% Process DFT structures for each symprec value
    dft_analysis_dict: dict[float, pd.DataFrame] = {}
    for symprec in symprec_values:
        symprec_str = f"symprec={symprec:.0e}".replace("e-0", "e-")

        # Always use full DFT analysis file, regardless of debug mode, ensuring
        # reference data available for all model structures regardless of debug mode
        # and sorting of material IDs
        dft_csv_path = (
            f"{ROOT}/data/wbm/dft-geo-opt-{symprec_str}-{moyo_version}.csv.gz"
        )

        if os.path.isfile(dft_csv_path):
            print(f"Loading DFT analysis from {dft_csv_path}")
            dft_analysis_dict[symprec] = pd.read_csv(dft_csv_path).set_index(Key.mat_id)
        else:
            dft_analysis_dict[symprec] = symmetry.get_sym_info_from_structs(
                dft_structs,
                pbar=dict(desc=f"Getting DFT symmetries {symprec=}"),
                symprec=symprec,
            )
            dft_analysis_dict[symprec].to_csv(dft_csv_path)

    # Create list of all model-symprec combinations
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
        # Wait for all tasks to complete
        for future in futures:
            future.result()  # This will raise any exceptions that occurred

    print(f"\nAll {len(tasks)} model-symprec combinations processed!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
