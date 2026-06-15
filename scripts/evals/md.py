"""Aggregate per-system MD metric files and write model-level metrics to YAML."""

import os
from glob import glob

import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.cli import cli_args
from matbench_discovery.md import default_md_reference_path, list_reference_systems
from matbench_discovery.metrics import md as md_metrics


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
        # only a dict carries pred_file/pred_file_url; a missing key (None) or a
        # 'not available' placeholder string must not flow into md_yaml.get(...) below
        md_yaml = model.metrics.get("md")
        md_yaml = md_yaml if isinstance(md_yaml, dict) else {}
        # collect per-system MD metric CSVs a parallel run wrote for a model, reading
        # each once and keeping single-system outputs (one row), excluding the
        # multi-system subset files that share this filename shape
        arch_dir = os.path.dirname(model.rel_path)
        pattern = f"{ROOT}/models/{arch_dir}/*md-nvt*/{model.name}-md-metrics-*.csv.gz"
        per_system_dfs = [
            df for path in sorted(glob(pattern)) if len(df := pd.read_csv(path)) == 1
        ]
        try:
            if per_system_dfs:  # parallel runs: concat per-system rows into one CSV
                df_md = md_metrics.combine_per_system_metrics(per_system_dfs)
                pred_file = f"models/{arch_dir}/{today}-{model.name}-md-metrics.csv.gz"
            else:  # fall back to an existing combined prediction file
                md_path = model.md_path  # getter may download the file
                if not isinstance(md_path, str) or not os.path.isfile(md_path):
                    print(f"Skipping {model.label}: no per-system CSVs or md_path")
                    n_skipped += 1
                    continue
                df_md = pd.read_csv(md_path)
                if "system" in df_md:  # index by system to match the per-system path
                    df_md = df_md.set_index("system")
                pred_file = md_yaml.get("pred_file") or md_path

            # don't persist model metrics from incomplete coverage (uneven system means
            # across models corrupt the leaderboard); guards both freshly combined
            # per-system files and a submitted combined CSV
            if expected is None:
                # canonical CFPMD-26 system names (authoritative coverage set);
                # triggers the reference download if not cached
                expected = set(list_reference_systems(default_md_reference_path()))
            if missing := expected - set(df_md.index):
                print(
                    f"Skipping {model.label}: {len(missing)}/{len(expected)} "
                    f"systems missing, e.g. {sorted(missing)[:3]}"
                )
                n_skipped += 1
                continue

            if per_system_dfs:  # persist the freshly combined per-system metrics
                df_md.to_csv(f"{ROOT}/{pred_file}")
                print(f"\n{model.label}: combined {len(per_system_dfs)} systems")

            metrics = md_metrics.calc_md_metrics(df_md)
            for key, value in metrics.items():
                shown = f"{value:.4f}" if isinstance(value, float) else value
                print(f"\t{key}={shown}")
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
    raise SystemExit(main())
