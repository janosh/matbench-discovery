"""Aggregate per-system MD metric files and write model-level metrics to YAML."""

import os
import sys
from glob import glob

import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.cli import cli_args
from matbench_discovery.enums import Model
from matbench_discovery.md import default_md_reference_paths
from matbench_discovery.metrics import md as md_metrics


def expected_systems() -> set[str]:
    """Canonical CFPMD-26 system directory names from the reference dataset, used to
    require complete coverage before writing model-level metrics to YAML. Triggers
    the reference download/extract if not cached (the authoritative system list).
    """
    ref_dir, _settings_csv = default_md_reference_paths()
    return {entry.name for entry in os.scandir(ref_dir) if entry.is_dir()}


def find_per_system_csvs(model: Model) -> list[str]:
    """Per-system MD metric CSVs a parallel run wrote for a model, across date dirs.

    Matches ``models/<arch>/<date>-md-nvt/<model>-md-metrics-<system>.csv.gz`` (the
    suffixed files single-system jobs produce), not the combined no-suffix CSV. The
    ``-md-metrics-`` separator keeps prefixes distinct (orb_v2 won't match
    orb_v2_mptrj). Multi-system subset outputs share this filename shape, so only
    keep one-row CSVs. Sorting puts later date dirs last so load_per_system_metrics'
    keep-last dedup lets a newer rerun override an earlier system result.
    """
    arch_dir = os.path.dirname(model.rel_path)
    pattern = f"{ROOT}/models/{arch_dir}/*md-nvt*/{model.name}-md-metrics-*.csv.gz"
    return [
        path
        for path in sorted(glob(pattern))
        if len(pd.read_csv(path, usecols=["system"])) == 1
    ]


def main() -> int:
    """Evaluate MD metrics and update model YAML files.

    Returns:
        Exit code: 0 if at least one model was evaluated, 1 otherwise.
    """
    models_to_evaluate = cli_args.models
    print(f"Evaluating MD metrics for {len(models_to_evaluate)} model(s)...")

    n_success = 0
    n_skipped = 0
    expected: set[str] | None = None  # CFPMD-26 system set, resolved lazily once

    for model in models_to_evaluate:
        md_yaml = model.metrics.get("md") or {}
        per_system_csvs = find_per_system_csvs(model)
        try:
            if per_system_csvs:  # parallel runs: concat per-system rows into one CSV
                df_md = md_metrics.load_per_system_metrics(per_system_csvs)
                # Don't persist model metrics from incomplete coverage: uneven means
                # across models corrupt the leaderboard.
                if expected is None:
                    expected = expected_systems()
                if missing := expected - set(df_md.index):
                    print(
                        f"Skipping {model.label}: {len(missing)}/{len(expected)} "
                        f"systems missing, e.g. {sorted(missing)[:3]}"
                    )
                    n_skipped += 1
                    continue
                arch_dir = os.path.dirname(model.rel_path)
                pred_file = f"models/{arch_dir}/{today}-{model.name}-md-metrics.csv.gz"
                df_md.to_csv(f"{ROOT}/{pred_file}")
                print(f"\n{model.label}: combined {len(per_system_csvs)} systems")
            else:  # fall back to an existing combined prediction file
                md_path = model.md_path  # getter may download the file
                if not isinstance(md_path, str) or not os.path.isfile(md_path):
                    print(f"Skipping {model.label}: no per-system CSVs or md_path")
                    n_skipped += 1
                    continue
                df_md = pd.read_csv(md_path)
                pred_file = md_yaml.get("pred_file") or md_path

            metrics = md_metrics.calc_md_metrics(df_md)
            for key, value in metrics.items():
                print(f"\t{key}={value:.4f}")
            md_metrics.write_metrics_to_yaml(
                model,
                metrics,
                pred_file_path=pred_file,
                pred_file_url=md_yaml.get("pred_file_url"),
            )
            print(f"\tUpdated {model.yaml_path}")
            n_success += 1
        except (ValueError, OSError, KeyError) as exc:
            print(f"\tError processing {model.label}: {exc}")
            n_skipped += 1

    if n_success == 0:
        print(f"\nNo models evaluated successfully ({n_skipped} skipped)")
        return 1
    print(f"\nSuccessfully evaluated {n_success} model(s), {n_skipped} skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
