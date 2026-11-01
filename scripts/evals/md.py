"""Aggregate per-system MD metric files and write model-level metrics to YAML."""

import os
from glob import escape, glob

import pandas as pd

from matbench_discovery import ROOT, today
from matbench_discovery.cli import cli_parser, parse_cli_args
from matbench_discovery.data import artifact_filename
from matbench_discovery.enums import Model
from matbench_discovery.md import default_md_reference_path, list_reference_systems
from matbench_discovery.metrics import md as md_metrics
from scripts.evals import evaluate_models

md_group = cli_parser.add_argument_group("md")
md_group.add_argument(
    "--md-run-dir",
    help="Combine MD metric CSVs from this directory instead of the declared file.",
)
md_group.add_argument(
    "--replace-author-reported",
    action="store_true",
    help="Let artifact-backed cost replace an author_reported YAML cost block.",
)


def resolve_metrics(
    model: Model, run_dir: str | None
) -> tuple[pd.DataFrame, str | None] | None:
    """Read the declared prediction or the explicitly selected run directory.

    Return a destination only for a fresh combination that the caller must persist.
    """
    pred_file = None
    if run_dir is None:
        # getter verifies/downloads the declared file and raises if it is missing
        if (md_path := model.md_path) is None:
            return None
        paths = [md_path]
    else:
        paths = sorted(glob(f"{escape(run_dir)}/*{model.name}-md-metrics-*.csv.gz"))
        if not paths:
            raise FileNotFoundError(f"No {model.name} MD metric CSVs in {run_dir!r}")
        pred_file = f"{model.artifact_dir}/{artifact_filename(today, 'md_metrics')}"
    # Retain duplicates so the coverage check rejects ambiguous reruns.
    return pd.concat(
        [
            pd.read_csv(path, index_col="system", float_precision="round_trip")
            for path in paths
        ]
    ), pred_file


def main(raw_args: list[str] | None = None) -> int:
    """Evaluate MD metrics with the shared sweep's failure and skip handling."""
    args = parse_cli_args(raw_args)
    expected: set[str] | None = None

    def evaluate_one(model: Model) -> str | None:
        nonlocal expected
        resolved = resolve_metrics(model, args.md_run_dir)
        if resolved is None:
            return "no declared MD prediction"
        df_md, new_pred_file = resolved

        # Resolve the canonical set only after finding an artifact.
        if expected is None:
            expected = set(list_reference_systems(default_md_reference_path()))
        if problems := md_metrics.coverage_problems(df_md.index, expected):
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
        md_metrics.write_metrics_to_yaml(
            model,
            metrics,
            pred_file_path=new_pred_file,
            replace_author_reported=args.replace_author_reported,
        )
        print(f"\tUpdated {model.yaml_path}")
        return None

    return evaluate_models("MD", args.models, evaluate_one)


if __name__ == "__main__":
    raise SystemExit(main())
