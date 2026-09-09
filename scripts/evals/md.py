"""Aggregate per-system MD metric files and write model-level metrics to YAML."""

import os
from glob import escape, glob

import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.cli import cli_args
from matbench_discovery.data import artifact_filename
from matbench_discovery.enums import Model
from matbench_discovery.md import default_md_reference_path, list_reference_systems
from matbench_discovery.metrics import md as md_metrics
from scripts.evals import evaluate_models


def resolve_metrics(
    model: Model, run_dir: str | None
) -> tuple[pd.DataFrame, str | None] | None:
    """Read the declared prediction or the explicitly selected run directory.

    Return a destination only for a fresh combination that the caller must persist.
    """
    pred_file = None
    if run_dir is None:
        md_path = model.md_path  # getter verifies/downloads the declared file
        if not md_path or not os.path.isfile(md_path):
            return None
        paths = [md_path]
    else:
        paths = sorted(glob(f"{escape(run_dir)}/*{model.name}-md-metrics-*.csv.gz"))
        if not paths:
            raise FileNotFoundError(f"No {model.name} MD metric CSVs in {run_dir!r}")
        model_dir = os.path.splitext(model.rel_path)[0]
        pred_file = f"models/{model_dir}/{artifact_filename(today, 'md_metrics')}"
    # Retain duplicates so the coverage check rejects ambiguous reruns.
    return pd.concat(
        [
            pd.read_csv(path, index_col="system", float_precision="round_trip")
            for path in paths
        ]
    ), pred_file


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
    """Evaluate MD metrics with the shared sweep's failure and skip handling."""
    expected: set[str] | None = None

    def evaluate_one(model: Model) -> str | None:
        nonlocal expected
        resolved = resolve_metrics(model, cli_args.md_run_dir)
        if resolved is None:
            return "no declared MD prediction"
        df_md, new_pred_file = resolved

        # Resolve the canonical set only after finding an artifact.
        if expected is None:
            expected = set(list_reference_systems(default_md_reference_path()))
        if problems := coverage_problems(df_md.index, expected):
            raise ValueError("; ".join(problems))

        metrics = md_metrics.calc_md_metrics(df_md)
        if new_pred_file:
            out_csv = f"{ROOT}/{new_pred_file}"
            os.makedirs(os.path.dirname(out_csv), exist_ok=True)
            df_md.to_csv(out_csv)
            print(f"\n{model.label}: combined {len(df_md)} systems")

        for key, value in metrics.items():
            shown = f"{value:.4f}" if isinstance(value, float) else value
            print(f"\t{key}={shown}")
        # Recomputing an existing file must retain its URL, size and checksum.
        md_metrics.write_metrics_to_yaml(model, metrics, pred_file_path=new_pred_file)
        print(f"\tUpdated {model.yaml_path}")
        return None

    return evaluate_models("MD", cli_args.models, evaluate_one)


if __name__ == "__main__":
    raise SystemExit(main())
