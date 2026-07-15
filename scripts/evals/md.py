"""Aggregate per-system MD metric files and write model-level metrics to YAML."""

import os
from glob import glob
from typing import Any

import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.cli import cli_args
from matbench_discovery.data import artifact_filename, file_ref_name, file_ref_url
from matbench_discovery.enums import Model
from matbench_discovery.md import default_md_reference_path, list_reference_systems
from matbench_discovery.metrics import md as md_metrics


def resolve_metrics(
    model: Model, md_yaml: dict[str, Any]
) -> tuple[pd.DataFrame, str, str | None, bool] | None:
    """Load a model's per-system MD metrics indexed by system, or None to skip.

    Prefers per-system CSVs from parallel single-system runs (combined into one frame);
    falls back to the model's submitted combined CSV (``md_path``). Returns
    ``(df_md, pred_file, pred_file_url, is_fresh_combine)``; ``is_fresh_combine`` flags
    the combined-from-per-system case, whose CSV the caller persists after the coverage
    check (and which carries no ``pred_file_url`` yet, so a stale one isn't reused).
    """
    arch_dir = os.path.dirname(model.rel_path)
    # one row each, excluding multi-system subset files that share this filename shape
    pattern = f"{ROOT}/models/{arch_dir}/*md-nvt*/*{model.name}-md-metrics-*.csv.gz"
    per_system_dfs = [
        df for path in sorted(glob(pattern)) if len(df := pd.read_csv(path)) == 1
    ]
    if per_system_dfs:  # parallel runs: concat per-system rows into one CSV
        df_md = md_metrics.combine_per_system_metrics(per_system_dfs)
        model_dir = os.path.splitext(model.rel_path)[0]
        pred_file = f"models/{model_dir}/{artifact_filename(today, 'md_metrics')}"
        return df_md, pred_file, None, True

    md_path = model.md_path  # getter may download the file
    if not md_path or not os.path.isfile(md_path):
        return None
    df_md = pd.read_csv(md_path)
    if "system" in df_md:  # index by system to match the per-system path
        df_md = df_md.set_index("system")
    pred_file = file_ref_name(md_yaml.get("pred_file")) or md_path
    return df_md, pred_file, file_ref_url(md_yaml.get("pred_file")), False


def coverage_problems(index: pd.Index, expected: set[str]) -> list[str]:
    """Reasons a per-system metric index isn't exactly the expected DynaMat v1.0 set
    (duplicate, missing or unexpected systems); an empty list means exact coverage.
    """
    present = set(index)
    problems = []
    if index.has_duplicates:
        dups = sorted(index[index.duplicated()].unique())
        problems.append(f"duplicate systems {dups}")
    if missing := expected - present:
        problems.append(f"{len(missing)} missing e.g. {sorted(missing)[:3]}")
    if extra := present - expected:
        problems.append(f"{len(extra)} unexpected e.g. {sorted(extra)[:3]}")
    return problems


def main() -> int:
    """Evaluate MD metrics and update model YAML files.

    Returns:
        Exit code: 0 if at least one model was evaluated, 1 otherwise.
    """
    models_to_evaluate = cli_args.models
    print(f"Evaluating MD metrics for {len(models_to_evaluate)} model(s)...")

    n_success = 0
    n_skipped = 0
    expected: set[str] | None = None  # DynaMat v1.0 system set, resolved lazily once

    for model in models_to_evaluate:
        md_yaml = model.metrics.get("md") or {}
        try:
            # resolve inside the try so a malformed CSV skips only this model
            resolved = resolve_metrics(model, md_yaml)
            if resolved is None:
                print(f"Skipping {model.label}: no per-system CSVs or md_path")
                n_skipped += 1
                continue
            df_md, pred_file, pred_file_url, is_fresh_combine = resolved

            # only persist metrics from exact DynaMat v1.0 coverage: missing systems
            # give uneven means across models, while duplicate or unexpected (e.g.
            # debug) rows silently skew them. Guards fresh and submitted CSVs alike.
            if expected is None:
                # canonical DynaMat v1.0 names (downloads reference if not cached)
                expected = set(list_reference_systems(default_md_reference_path()))
            if problems := coverage_problems(df_md.index, expected):
                print(f"Skipping {model.label}: " + "; ".join(problems))
                n_skipped += 1
                continue

            if is_fresh_combine:  # persist the freshly combined per-system metrics
                out_csv = f"{ROOT}/{pred_file}"
                os.makedirs(os.path.dirname(out_csv), exist_ok=True)
                df_md.to_csv(out_csv)
                print(f"\n{model.label}: combined {len(df_md)} systems")

            metrics = md_metrics.calc_md_metrics(df_md)
            for key, value in metrics.items():
                shown = f"{value:.4f}" if isinstance(value, float) else value
                print(f"\t{key}={shown}")
            md_metrics.write_metrics_to_yaml(
                model, metrics, pred_file_path=pred_file, pred_file_url=pred_file_url
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
