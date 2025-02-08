"""Run this script to add/update geometry optimization analysis for new models to
individual CSV files (one per symprec value) in a model's directory.

Output files will have the same name as the input file containing the model's relaxed
structures, but with the symprec value appended to the filename.

Example usage:
    python scripts/analyze_geo_opt.py --models mace_mp_0 m3gnet --symprec 1e-2 1e-5
    python scripts/analyze_geo_opt.py --debug 10  # only analyze first 10 structures
    python scripts/analyze_geo_opt.py --workers 4  # use 4 CPU cores
"""

# %%
import argparse
import importlib
import importlib.metadata
import itertools
import multiprocessing as mp
import os
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Final

import pandas as pd
from pymatgen.core import Structure
from pymatviz.enums import Key

from matbench_discovery import ROOT
from matbench_discovery.enums import DataFiles, Model
from matbench_discovery.metrics import geo_opt
from matbench_discovery.models import MODEL_METADATA
from matbench_discovery.structure import analyze_symmetry, pred_vs_ref_struct_symmetry


def analyze_model_symprec(
    model_name: str,
    symprec: float,
    df_dft_analysis: pd.DataFrame,
    dft_structs: dict[str, Structure],
    debug_mode: int = 0,
    moyo_version: str = "moyo=0.3.1",
    pbar_pos: int = 0,  # tqdm progress bar position
) -> None:
    """Analyze a single model for a single symprec value."""
    model = Model[model_name]
    model_metadata = MODEL_METADATA[model.label]

    geo_opt_metrics: dict[str, Any] = model_metadata.get("metrics", {}).get(
        "geo_opt", {}
    )

    # skip models that don't support geometry optimization
    if geo_opt_metrics in ("not applicable", "not available"):
        print(f"⚠️ {model.label} does not support geometry optimization")
        return

    if not model.geo_opt_path:
        print(f"⚠️ {model.label} has no relaxed structures file")
        return

    if not os.path.isfile(ml_relaxed_structs_path := model.geo_opt_path):
        print(
            f"⚠️ {model.label}-relaxed structures not found, expected "
            f"at {ml_relaxed_structs_path}"
        )
        return

    # Load model structures
    try:
        if ml_relaxed_structs_path.endswith((".json", ".json.gz", ".json.xz")):
            df_ml_structs = pd.read_json(ml_relaxed_structs_path)
        else:
            raise ValueError(
                "Relaxed structure analysis currently only supports pymatgen JSON, "
                f"got {ml_relaxed_structs_path}"
            )
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
        return

    # Convert structures
    model_structs = {
        mat_id: Structure.from_dict(struct_dict)
        for mat_id, struct_dict in df_ml_structs[struct_col].items()
    }

    symprec_str = f"symprec={symprec:.0e}".replace("e-0", "e-")
    geo_opt_filename = model.geo_opt_path.removesuffix(".json.gz")
    geo_opt_csv_path = f"{geo_opt_filename}-{symprec_str}-{moyo_version}.csv.gz"

    if os.path.isfile(geo_opt_csv_path):
        print(f"{model.label} already analyzed at {geo_opt_csv_path}")
        return

    # Analyze symmetry for current symprec
    pbar_desc = f"Process {pbar_pos}: Analyzing {model.label} for {symprec=}"
    df_model_analysis = analyze_symmetry(
        model_structs,
        pbar=dict(desc=pbar_desc, position=pbar_pos, leave=True),
        symprec=symprec,
    )

    # Compare with DFT reference
    pbar_desc = f"Process {pbar_pos}:Comparing DFT vs {model.label} for {symprec=}"
    # break here
    df_ml_geo_analysis = pred_vs_ref_struct_symmetry(
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
    print(f"\nCalculated metrics: {df_metrics}")  # Debug print
    geo_opt.write_geo_opt_metrics_to_yaml(df_metrics, model, symprec)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--models",
        nargs="+",
        default=["all"],
        help="Model names to analyze. Use 'all' to analyze all models.",
    )
    parser.add_argument(
        "--symprec",
        nargs="+",
        type=float,
        default=[1e-2, 1e-5],
        help="Symmetry precision values to analyze.",
    )
    parser.add_argument(
        "--debug",
        type=int,
        default=0,
        help="If > 0, only analyze this many structures.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, mp.cpu_count() - 1),
        help="Number of processes to use. Defaults to number of model-symprec combos.",
    )
    args = parser.parse_args()

    # set to > 0 to activate debug mode, only that many structures will be analyzed
    debug_mode: Final[int] = args.debug
    # List of symprec values to analyze
    symprec_values: Final[Sequence[float]] = args.symprec

    # Get list of models to analyze
    if "all" in args.models:
        model_names = [model.name for model in Model]
    else:
        model_names = args.models
        # Validate model names
        valid_models = {model.name for model in Model}
        for model_name in model_names:
            if model_name not in valid_models:
                raise ValueError(f"Invalid {model_name=}, valid models: {valid_models}")

    moyo_version = f"moyo={importlib.metadata.version('moyopy')}"

    # %% Load WBM reference structures (this takes a while)
    print("Loading WBM reference structures...")
    wbm_cse_path = DataFiles.wbm_computed_structure_entries.path
    df_wbm_structs: pd.DataFrame = pd.read_json(wbm_cse_path).set_index(Key.mat_id)

    if debug_mode:
        df_wbm_structs = df_wbm_structs.head(debug_mode)

    dft_structs: dict[str, Structure] = {
        mat_id: Structure.from_dict(cse[Key.structure])
        for mat_id, cse in df_wbm_structs[Key.computed_structure_entry].items()
    }

    # %% Process DFT structures for each symprec value
    dft_analysis_dict: dict[float, pd.DataFrame] = {}
    for symprec in symprec_values:
        symprec_str = f"symprec={symprec:.0e}".replace("e-0", "e-")
        dft_csv_path = (
            f"{ROOT}/data/wbm/dft-geo-opt-{symprec_str}-{moyo_version}.csv.gz"
        )

        if os.path.isfile(dft_csv_path):
            dft_analysis_dict[symprec] = pd.read_csv(dft_csv_path).set_index(Key.mat_id)
        else:
            dft_analysis_dict[symprec] = analyze_symmetry(
                dft_structs,
                pbar=dict(desc=f"Getting DFT symmetries {symprec=}"),
                symprec=symprec,
            )
            dft_analysis_dict[symprec].to_csv(dft_csv_path)

    # Create list of all model-symprec combinations
    tasks = list(itertools.product(model_names, symprec_values))
    n_workers = min(len(tasks), args.workers)

    # %% Process model-symprec combinations in parallel
    print(
        f"\nAnalyzing {len(tasks)} model-symprec combos using {n_workers} processes..."
    )

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(
                analyze_model_symprec,
                model_name,
                symprec,
                dft_analysis_dict[symprec],
                dft_structs,
                debug_mode,
                moyo_version,
                pbar_pos=idx,  # assign unique position to each task's progress bar
            )
            for idx, (model_name, symprec) in enumerate(tasks)
        ]
        # Wait for all tasks to complete
        for future in futures:
            future.result()  # This will raise any exceptions that occurred

    print(f"\nAll {len(tasks)} model-symprec combinations processed!")
