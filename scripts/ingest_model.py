"""Model-submission ingestion pipeline.

The checklist reads the model YAML through the Model enum's metadata API instead of
grep/sed, making it testable (tests/test_ingest_model.py) and robust to YAML
formatting.

Usage:
    uv run scripts/ingest_model.py <model>            # check+evals+figs+payloads
    uv run scripts/ingest_model.py <model> --archive  # + figshare/release upload
    uv run scripts/ingest_model.py --payloads-only    # full payload regen
    uv run scripts/ingest_model.py <model> --payloads-only  # merge single model
"""

import argparse
import glob
import os
import shlex
import subprocess
from collections import Counter
from collections.abc import Sequence
from functools import partialmethod
from pathlib import Path

from matbench_discovery import ROOT
from matbench_discovery.calculators import CALCULATORS
from matbench_discovery.discovery import ARCHIVED_DISCOVERY_MODELS, RelaxationSettings
from matbench_discovery.enums import Model

PASS, FAIL, SKIP = "✓", "✗", "○"
PAYLOAD_FLAGS = ("--auto-download",)
# each entry is the `uv run` argument string for one payload script; kappa needs the
# phonons extra (phono3py/phonopy) when computing conductivity diagnostics
PAYLOAD_SCRIPTS = (
    "scripts/model_figs/roc_curves_models.py",
    "scripts/model_figs/hull_dist_box_plot.py",
    "scripts/model_figs/cumulative_metrics.py",
    "scripts/model_figs/rolling_hull_dist_mae_models.py",
    "scripts/model_figs/tiles_hist_classified_stable_models.py",
    "scripts/model_figs/tmi-page-figures.py",
    "scripts/model_figs/single_model_per_element_errors.py",
    "--extra phonons scripts/model_figs/kappa_103_analysis.py",
    "scripts/evals/geo_opt.py",
)
PARITY_ASSET_DIRS = (
    "site/static/energy-parity/assets",
    "site/static/kappa-parity/assets",
)
# (label, skipped for energy-only models, fatal on failure, `uv run` args).
# kappa needs the phonons extra (phono3py/phonopy); geo-opt the symmetry extra (moyopy)
EVAL_STEPS = (
    ("Discovery metrics", False, True, "scripts/evals/discovery.py"),
    ("Kappa metrics", True, True, "--extra phonons scripts/evals/kappa.py"),
    ("Geo-opt analysis", True, True, "--extra symmetry scripts/analyze_geo_opt.py"),
    ("Diatomic metrics", True, False, "scripts/evals/diatomic_metrics.py"),
)
FIG_STEPS = (
    (
        "Energy parity assets",
        False,
        True,
        "site/scripts/generate-energy-parity-assets.py --skip-structures",
    ),
    (
        "Kappa parity assets",
        True,
        True,
        "site/scripts/generate-kappa-parity-assets.py",
    ),
)


class Checklist:
    """Collects pass/fail/skip results and renders a summary."""

    def __init__(self) -> None:
        self.results: list[tuple[str, str]] = []

    def record(self, status: str, msg: str) -> None:
        """Append a (status, message) result and print it to the console."""
        self.results.append((status, msg))
        suffix = " (skipped/optional)" if status == SKIP else ""
        print(f"  {status} {msg}{suffix}")

    ok = partialmethod(record, PASS)
    fail = partialmethod(record, FAIL)
    skip = partialmethod(record, SKIP)

    @property
    def n_failed(self) -> int:
        """Number of recorded results with FAIL status."""
        return sum(status == FAIL for status, _ in self.results)

    def summary(self) -> str:
        """Render the passed/failed/skipped counts as a multi-line summary string."""
        counts = Counter(status for status, _ in self.results)
        return (
            f"  {PASS} Passed:  {counts[PASS]}\n"
            f"  {FAIL} Failed:  {counts[FAIL]}\n"
            f"  {SKIP} Skipped: {counts[SKIP]}"
        )


def banner(text: str) -> None:
    print(f"\n{'━' * 78}\n{text}\n{'━' * 78}")


def run_cmd(*cmd: str) -> bool:
    """Run a command, streaming output; return True on zero exit."""
    print(f">> {' '.join(cmd)}")
    return subprocess.run(cmd, check=False).returncode == 0


def uv_run_args(args: str) -> tuple[str, ...]:
    """Build ``uv run`` args, resolving leading ``--extra NAME`` flags into an
    editable project requirement (e.g. '--extra phonons x.py' -> --with-editable
    .[phonons] x.py).
    """
    tokens = shlex.split(args)
    extras: list[str] = []
    while tokens[:1] == ["--extra"]:
        if len(tokens) < 2 or not tokens[1]:
            raise ValueError("uv run --extra requires a non-empty value")
        extras.append(tokens[1])
        tokens = tokens[2:]
    project_req = f".[{','.join(extras)}]" if extras else "."
    return ("uv", "run", "--with-editable", project_req, *tokens)


def task_metrics(model: Model, task: str) -> dict | str | None:
    """A task's raw metrics from the model YAML: dict if present, str for declared
    opt-outs like 'not available', None if missing.
    """
    return model.metadata.get("metrics", {}).get(task)


def check_submission(model: Model, checks: Checklist) -> bool:
    """Validate PR checklist requirements via the model's parsed YAML metadata.

    Returns True when the model is energy-only (targets=E), which downstream steps
    use to skip force-based tasks (geo-opt, phonons, diatomics).
    """
    banner("STEP 1: Checking PR checklist requirements")
    yaml_path = Path(model.yaml_path)
    if yaml_path.is_file():
        checks.ok(f"Model YAML file exists: {yaml_path}")
    else:
        checks.fail(f"Model YAML file not found at: {yaml_path}")
    checks.ok(f"Model '{model.name}' is registered in the Model enum")

    # E = energy only (no forces -> no relaxation, phonons or diatomics)
    energy_only = model.metadata.get("targets") == "E"
    calc_spec = CALCULATORS.get(model.name)

    discovery = task_metrics(model, "discovery")
    if isinstance(discovery, dict) and discovery.get("pred_file_url"):
        checks.ok("Prediction file URL found in YAML (discovery.pred_file_url)")
    else:
        checks.fail(f"No discovery.pred_file_url found in {yaml_path}")

    # (task, label, optional). phonons pred file nests under kappa_103.
    force_tasks = (
        ("geo_opt", "Relaxed structures file (geo_opt.pred_file)", False),
        ("phonons", "Phonon predictions (phonons.kappa_103.pred_file)", False),
        ("diatomics", "Diatomics predictions (diatomics.pred_file)", True),
    )
    for task, label, optional in force_tasks:
        if energy_only:
            checks.skip(f"{label} check skipped (targets=E, no forces)")
            continue
        metrics = task_metrics(model, task)
        if isinstance(metrics, str):  # declared opt-out, e.g. 'not available'
            checks.skip(f"{label} declared {metrics!r}")
            continue
        if task == "phonons" and isinstance(metrics, dict):
            metrics = metrics.get("kappa_103")
        if isinstance(metrics, dict) and metrics.get("pred_file"):
            checks.ok(f"{label} found in YAML")
        elif optional:
            checks.skip(f"{label} not found in YAML")
        else:
            checks.fail(f"{label} not found in {yaml_path}")

    for task in ("discovery", "kappa", "diatomics"):
        if energy_only and task != "discovery":
            checks.skip(f"{task} test script check skipped (targets=E, no forces)")
            continue
        if task == "discovery" and (
            reason := ARCHIVED_DISCOVERY_MODELS.get(model.name)
        ):
            checks.skip(f"discovery is archived: {reason}")
            continue
        if task == "discovery" and calc_spec is None:
            checks.fail(
                "discovery model is not registered with the shared runner: add an "
                "ASE calculator to matbench_discovery/calculators.py, or discuss "
                "non-ASE inference with the maintainers first (see contributing.md)"
            )
            continue

        scripts = sorted(yaml_path.parent.glob(f"test_*_{task}.py"))
        if task != "discovery" and scripts:
            checks.ok(f"{task} test script found: {scripts[0].name}")
        elif task in ("discovery", "diatomics") and calc_spec is not None:
            runner = f"{ROOT}/models/run_{task}.py"
            try:
                if not os.path.isfile(runner):
                    raise ValueError(f"shared runner not found: {runner}")
                calc_spec.uv_run_cmd(runner, "--model", model.name, "--dry-run")
                if task == "discovery":
                    RelaxationSettings.from_model(model.name)
            except (TypeError, ValueError) as exc:
                checks.fail(f"Invalid shared {task} runner configuration: {exc}")
            else:
                checks.ok(f"{task} uses shared runner: models/run_{task}.py")
        elif task == "diatomics":
            checks.skip(f"{task} test script not found (test_*_{task}.py)")
        else:
            checks.fail(
                f"{task} test script not found (test_*_{task}.py) in {yaml_path.parent}"
            )

    return energy_only


def run_model_steps(
    title: str,
    steps: Sequence[tuple[str, bool, bool, str]],
    model: Model,
    checks: Checklist,
    *,
    energy_only: bool,
    extra_args: Sequence[str] = (),
) -> None:
    """Run a table of per-model commands, recording pass/fail/skip for each."""
    banner(title)
    for label, needs_forces, fatal, args in steps:
        if needs_forces and energy_only:
            checks.skip(f"{label} skipped (targets=E, no forces)")
        elif run_cmd(*uv_run_args(args), "--models", model.name, *extra_args):
            checks.ok(f"{label} completed")
        elif fatal:
            checks.fail(f"{label} failed")
        else:
            checks.skip(f"{label} skipped (no data or failed)")


def run_archive(model: Model, checks: Checklist) -> None:
    """Archive the model's prediction files to the project's figshare articles (one
    per prediction task) for longevity - rewrites the YAML's *_url keys - and publish
    new parity assets to the GitHub release the site build downloads from.
    """
    banner("STEP 4: Archiving prediction files + publishing parity assets")
    if run_cmd(
        *uv_run_args("scripts/upload_model_preds_to_figshare.py"),
        *("--models", model.name, "--no-interactive"),
    ):
        checks.ok("Prediction files archived to project figshare articles")
    else:
        checks.fail("Figshare archival failed")

    for assets_dir in PARITY_ASSET_DIRS:
        if not (assets := glob.glob(f"{assets_dir}/*.json.gz")):
            continue
        if run_cmd("gh", "release", "upload", "v1.0.0", *assets, "--clobber"):
            checks.ok(f"{len(assets)} parity assets from {assets_dir} published")
        else:
            checks.fail(f"Publishing parity assets from {assets_dir} failed")


def run_payload_refresh(checks: Checklist, model: Model | None = None) -> None:
    """Refresh site/src/figs/*.jsonl plus route-local JSONL payloads, then test.

    With a model given, payload scripts run with --models <model> and splice only that
    model's freshly computed entries into the committed payloads (no other model's
    prediction files needed, see figs.write_site_payload); without one, payloads are
    regenerated from the full active roster.
    """
    suffix = f" for {model.name}" if model else ""
    banner(f"Refreshing multi-model site figure payloads{suffix}")
    model_args = ("--models", model.name) if model else ()
    for script in PAYLOAD_SCRIPTS:
        if not run_cmd(*uv_run_args(script), *PAYLOAD_FLAGS, *model_args):
            checks.fail(f"{script} failed")
            return
    if run_cmd(*uv_run_args("--with pytest pytest tests/test_fig_payloads.py -q")):
        checks.ok("Multi-model payloads regenerated + shape tests passed")
    else:
        checks.fail("Payload shape tests failed after regeneration")


def map_yaml_paths(paths: Sequence[str]) -> list[str]:
    """Map model YAML paths (as in a PR diff) to Model enum member names.

    Raises SystemExit for paths without an enum entry (used by CI to fail closed on
    submissions that forgot to register their model).
    """
    names: list[str] = []
    for path in paths:
        model = next((mdl for mdl in Model if mdl.yaml_path.endswith(path)), None)
        if model is None:
            raise SystemExit(f"No Model enum entry maps to {path} - add it to enums.py")
        names.append(model.name)
    return list(dict.fromkeys(names))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", nargs="?", help="Model enum name (dashes ok)")
    bool_flags = {
        "--overwrite": "Overwrite existing eval outputs",
        "--archive": "Also archive pred files to figshare + publish parity assets "
        "(needs FIGSHARE_TOKEN and gh auth)",
        "--payloads-only": "Only refresh the multi-model site figure payloads "
        "(merging just <model>'s entries if a model is given)",
    }
    for flag, help_msg in bool_flags.items():
        parser.add_argument(flag, action="store_true", help=help_msg)
    parser.add_argument(
        "--map-yaml-paths",
        nargs="+",
        metavar="PATH",
        help="Print enum names for model YAML paths (used by CI model detection)",
    )
    args = parser.parse_args(argv)

    if args.map_yaml_paths:
        print(" ".join(map_yaml_paths(args.map_yaml_paths)))
        return 0

    checks = Checklist()
    model: Model | None = None
    if args.model:
        try:  # Model(...) invokes _missing_ to normalize dashes/casing
            model = Model(args.model)
        except ValueError:
            parser.error(
                f"{args.model!r} not in Model enum - add it to "
                "matbench_discovery/enums.py"
            )
    if args.payloads_only:
        run_payload_refresh(checks, model=model)
    else:
        if model is None:
            parser.error("model is required unless --payloads-only")
        if args.archive and not os.getenv("FIGSHARE_TOKEN"):
            parser.error("FIGSHARE_TOKEN must be set for --archive")

        banner(f"Ingesting model submission: {model.name}")
        energy_only = check_submission(model, checks)
        overwrite = ("--overwrite",) if args.overwrite else ()
        run_model_steps(
            "STEP 2: Running evaluation scripts", EVAL_STEPS, model, checks,
            energy_only=energy_only, extra_args=("--auto-download", *overwrite),
        )  # fmt: skip
        run_model_steps(
            "STEP 3: Generating per-model figures", FIG_STEPS, model, checks,
            energy_only=energy_only,
        )  # fmt: skip
        if args.archive:
            run_archive(model, checks)
        run_payload_refresh(checks, model=model)

    banner("SUMMARY")
    print(checks.summary())
    if checks.n_failed:
        print(
            "\n⚠️  Some checks failed. Please review the PR checklist:\n"
            "   https://github.com/janosh/matbench-discovery/blob/main/.github/pull_request_template.md"
        )
        return 1
    print("\n✅ All required checks passed!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
