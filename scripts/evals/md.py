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
from scripts.evals import evaluate_models


def resolve_metrics(
    model: Model, md_yaml: dict[str, Any]
) -> tuple[pd.DataFrame, str, str | None, bool] | None:
    """Load per-system metrics when available, else the submitted combined CSV.

    The final return flag marks a fresh combination that the caller must persist.
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
    """Return duplicate, missing, and unexpected system coverage problems."""
    present = set(index)
    problems = []
    if index.has_duplicates:
        duplicates = sorted(index[index.duplicated()].unique())
        problems.append(f"duplicate systems {duplicates}")
    if missing := expected - present:
        problems.append(f"{len(missing)} missing e.g. {sorted(missing)[:3]}")
    if extra := present - expected:
        problems.append(f"{len(extra)} unexpected e.g. {sorted(extra)[:3]}")
    return problems


def main() -> int:
    """Evaluate MD metrics, returning 0 if any model succeeds and 1 otherwise."""
    expected: set[str] | None = None

    def evaluate_one(model: Model) -> str | None:
        nonlocal expected
        resolved = resolve_metrics(model, model.metrics.get("md") or {})
        if resolved is None:
            return "no per-system CSVs or md_path"
        df_md, pred_file, pred_file_url, is_fresh_combine = resolved

        # Resolve the canonical set only after finding an artifact.
        if expected is None:
            expected = set(list_reference_systems(default_md_reference_path()))
        if problems := coverage_problems(df_md.index, expected):
            return "; ".join(problems)

        if is_fresh_combine:
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
        return None

    return evaluate_models("MD", cli_args.models, evaluate_one)


if __name__ == "__main__":
    raise SystemExit(main())
