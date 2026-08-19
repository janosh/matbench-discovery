"""Model-submission ingestion pipeline.

The checklist reads the model YAML through the Model enum's metadata API instead of
grep/sed, making it testable (tests/scripts/test_ingest_model.py) and robust to YAML
formatting.

Usage:
    uv run scripts/ingest_model.py <model>            # check+evals+figs+payloads
    uv run scripts/ingest_model.py <model> --archive  # + figshare/release upload
    uv run scripts/ingest_model.py --payloads-only --full-roster
    uv run scripts/ingest_model.py <models...> --payloads-only  # merge models
"""

import argparse
import json
import os
import shlex
import shutil
import subprocess
from collections import Counter
from collections.abc import Sequence
from functools import partialmethod

import yaml

from matbench_discovery import ROOT
from matbench_discovery.calculators import CALCULATORS
from matbench_discovery.data import (
    file_ref_url,
    iter_file_refs,
    parse_artifact_filename,
    task_coverage,
)
from matbench_discovery.discovery import ARCHIVED_DISCOVERY_MODELS, RelaxationSettings
from matbench_discovery.enums import Model
from matbench_discovery.phonons.pipeline import KappaSettings

PASS, FAIL, SKIP = "✓", "✗", "○"
# task_coverage statuses that mean the task has results worth checking artifacts for
COVERED_STATUSES = frozenset({"complete", "partial"})
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


def check_submission(
    model: Model,
    checks: Checklist,
    *,
    validate_runner: bool = True,
) -> bool:
    """Validate the checklist and return whether the model predicts energy only.

    ``validate_runner=False`` limits trusted CI to submitted metadata and artifacts.
    """
    banner("STEP 1: Checking PR checklist requirements")
    yaml_path = model.yaml_path
    if os.path.isfile(yaml_path):
        checks.ok(f"Model YAML file exists: {yaml_path}")
    else:
        checks.fail(f"Model YAML file not found at: {yaml_path}")
    checks.ok(f"Model '{model.name}' is present in the generated Model enum")
    metadata = model.metadata
    metrics_by_task = model.metrics
    expected_artifact_dir = f"models/{model.family}/{metadata['model_key']}"
    try:
        artifact_refs = list(iter_file_refs(metrics_by_task, ("metrics",)))
    except ValueError as exc:
        checks.fail(str(exc))
        artifact_refs = []
    for key_path, artifact_path in artifact_refs:
        try:
            parse_artifact_filename(os.path.basename(artifact_path))
            if os.path.dirname(artifact_path) != expected_artifact_dir:
                raise ValueError("artifact is not directly model-owned")
        except ValueError as exc:
            checks.fail(f"Invalid {'.'.join(key_path)}: {exc}")

    # E = energy only (no forces -> no relaxation, phonons or diatomics)
    energy_only = metadata.get("targets") == "E"
    discovery = metrics_by_task.get("discovery")
    if isinstance(discovery, dict) and file_ref_url(discovery.get("pred_file")):
        checks.ok("Prediction file URL found in YAML (discovery.pred_file.url)")
    else:
        checks.fail(f"No discovery.pred_file.url found in {yaml_path}")

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
        status, reason = task_coverage(metadata, task)
        if status not in COVERED_STATUSES:
            suffix = f": {reason}" if reason else ""
            checks.skip(f"{label} declared {status!r}{suffix}")
            continue
        task_metrics = metrics_by_task.get(task)
        if task == "phonons" and isinstance(task_metrics, dict):
            task_metrics = task_metrics.get("kappa_103")
        if isinstance(task_metrics, dict) and task_metrics.get("pred_file"):
            checks.ok(f"{label} found in YAML")
        else:
            record = checks.skip if optional else checks.fail
            location = "YAML" if optional else yaml_path
            record(f"{label} not found in {location}")

    if not validate_runner:
        return energy_only

    calc_spec = CALCULATORS.get(model.name)
    phonons_available = task_coverage(metadata, "phonons")[0] in COVERED_STATUSES
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

        if task == "kappa" and not phonons_available:
            checks.skip("kappa shared runner check skipped (phonons unavailable)")
            continue

        if calc_spec is None:
            record = checks.fail if task == "kappa" else checks.skip
            record(
                f"{task} shared runner unsupported: no registered calculator for "
                f"{model.name}"
            )
            continue
        runner = f"{ROOT}/models/run_{task}.py"
        try:
            if not os.path.isfile(runner):
                raise ValueError(f"shared runner not found: {runner}")
            if not (task == "kappa" and calc_spec.requires_checkpoint):
                command = calc_spec.uv_run_cmd(
                    runner, "--model", model.name, "--dry-run"
                )
                subprocess.run(command, check=True)
            if task == "discovery":
                RelaxationSettings.from_model(model.name)
            elif task == "kappa":
                KappaSettings.from_model(model.name)
        except (
            subprocess.CalledProcessError,
            OSError,
            TypeError,
            ValueError,
        ) as exc:
            checks.fail(f"Invalid shared {task} runner configuration: {exc}")
        else:
            checks.ok(f"{task} uses shared runner: models/run_{task}.py")

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
    """Run model-specific commands, recording each result."""
    banner(title)
    model_args = ("--models", model.name)
    for label, needs_forces, fatal, args in steps:
        if needs_forces and energy_only:
            checks.skip(f"{label} skipped (targets=E, no forces)")
        elif run_cmd(*uv_run_args(args), *model_args, *extra_args):
            checks.ok(f"{label} completed")
        elif fatal:
            checks.fail(f"{label} failed")
        else:
            checks.skip(f"{label} skipped (no data or failed)")


def run_archive(model: Model, checks: Checklist) -> None:
    """Archive model predictions to project Figshare articles for longevity."""
    banner("STEP 4: Archiving prediction files")
    if run_cmd(
        *uv_run_args("scripts/upload_model_preds_to_figshare.py"),
        *("--models", model.name, "--no-interactive"),
    ):
        checks.ok("Prediction files archived to project figshare articles")
    else:
        checks.fail("Figshare archival failed")


def release_asset_digests() -> dict[str, str | None]:
    """Return each shared GitHub release asset's SHA-256 digest."""
    gh_path = shutil.which("gh")
    if gh_path is None:
        raise FileNotFoundError("GitHub CLI executable 'gh' not found")
    output = subprocess.check_output(
        (gh_path, "release", "view", "v1.0.0", "--json", "assets"),
        text=True,
    )
    return {
        asset["name"]: asset.get("digest") for asset in json.loads(output)["assets"]
    }


def publish_parity_assets(checks: Checklist) -> None:
    """Publish new parity assets without overwriting immutable release objects."""
    banner("Publishing parity assets")
    try:
        published = release_asset_digests()
    except (
        OSError,
        subprocess.CalledProcessError,
        json.JSONDecodeError,
        KeyError,
    ) as exc:
        checks.fail(f"Reading published parity assets failed: {exc}")
        return
    for parity_type in ("energy", "kappa"):
        assets_dir = f"site/static/{parity_type}-parity/assets"
        manifest_path = f"site/src/lib/parity/{parity_type}-parity-manifest.json"
        with open(manifest_path, encoding="utf-8") as file:
            manifest = json.load(file)
        entries = [
            manifest["base"],
            *manifest["model_assets"].values(),
            *manifest.get("structure_bundles", ()),
        ]
        pending: list[str] = []
        conflicts: list[str] = []
        for entry in entries:
            name = entry["asset"]
            asset_path = f"{assets_dir}/{name}"
            if name not in published and os.path.isfile(asset_path):
                pending.append(asset_path)
            elif published.get(name) != f"sha256:{entry['sha256']}":
                conflicts.append(name)
        if conflicts:
            checks.fail(
                f"{parity_type} parity assets missing locally or conflict with "
                f"immutable release objects: {', '.join(conflicts)}"
            )
        elif not pending:
            checks.ok(
                f"All {len(entries)} {parity_type} parity assets already published"
            )
        elif run_cmd("gh", "release", "upload", "v1.0.0", *sorted(pending)):
            checks.ok(f"{len(pending)} parity assets from {assets_dir} published")
        else:
            checks.fail(f"Publishing parity assets from {assets_dir} failed")


def run_payload_refresh(
    checks: Checklist,
    models: Sequence[Model] = (),
) -> None:
    """Refresh site/src/figs/*.jsonl plus route-local JSONL payloads, then test.

    With models given, payload scripts splice only those freshly computed entries into
    the committed payloads. Without models, they regenerate the full active roster.
    """
    suffix = f" for {', '.join(model.name for model in models)}" if models else ""
    banner(f"Refreshing multi-model site figure payloads{suffix}")
    payload_args = (
        ("--models", *(model.name for model in models))
        if models
        else ("--full-roster",)
    )
    for script in PAYLOAD_SCRIPTS:
        if not run_cmd(*uv_run_args(script), "--auto-download", *payload_args):
            checks.fail(f"{script} failed")
            return
    if not run_cmd(*uv_run_args("scripts/summarize_payload_changes.py")):
        checks.fail("Payload change report failed")
        return
    if run_cmd(*uv_run_args("--with pytest pytest tests/site/test_fig_payloads.py -q")):
        checks.ok("Multi-model payloads regenerated + shape tests passed")
    else:
        checks.fail("Payload shape tests failed after regeneration")


def read_head_yaml(path: str) -> dict[str, object] | None:
    """Read one model YAML from HEAD, returning None when it is newly added."""
    command = ["git", "show", f"HEAD:{path}"]
    result = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        stderr = result.stderr.strip()
        missing_path = f"path '{path}'"
        if missing_path in stderr and stderr.endswith(
            ("does not exist in 'HEAD'", "not in 'HEAD'")
        ):
            return None
        result.check_returncode()
    metadata = yaml.safe_load(result.stdout)
    if not isinstance(metadata, dict):
        raise TypeError(f"HEAD:{path} does not contain model metadata")
    return metadata


def classify_yaml_changes(
    changed_paths: Sequence[str], removed_paths: Sequence[str]
) -> dict[str, object]:
    """Classify model-YAML changes into targeted, full-roster, or no refresh."""
    new_by_key: dict[str, dict[str, object]] = {}
    old_by_key: dict[str, dict[str, object]] = {}
    for path in changed_paths:
        with open(f"{ROOT}/{path}", encoding="utf-8") as file:
            metadata = yaml.safe_load(file)
        if not isinstance(metadata, dict) or not isinstance(
            model_key := metadata.get("model_key"), str
        ):
            raise TypeError(f"{path} does not contain a string model_key")
        new_by_key[model_key] = metadata
        if (old_metadata := read_head_yaml(path)) is not None:
            old_key = old_metadata.get("model_key")
            if old_key != model_key:
                raise ValueError(
                    f"Automated ingestion cannot change model_key in {path}: "
                    f"{old_key!r} -> {model_key!r}; model keys are immutable"
                )
            old_by_key[model_key] = old_metadata
    for path in removed_paths:
        if (old_metadata := read_head_yaml(path)) is None:
            raise ValueError(f"Removed model YAML not found in HEAD: {path}")
        old_key = old_metadata.get("model_key")
        if not isinstance(old_key, str):
            raise TypeError(f"HEAD:{path} does not contain a string model_key")
        old_by_key[old_key] = old_metadata

    old_keys = set(old_by_key)
    removed_keys = old_keys - new_by_key.keys()
    lifecycle_changed = any(
        old_by_key[key].get("lifecycle") != new_by_key[key].get("lifecycle")
        for key in old_keys & new_by_key.keys()
    )
    return {
        "full_roster": bool(removed_keys or lifecycle_changed),
        "targeted_models": [
            Model.from_ref(key).name
            for key, metadata in sorted(new_by_key.items())
            if old_by_key.get(key) != metadata and metadata.get("lifecycle") == "active"
        ],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "models",
        nargs="*",
        type=Model.from_ref,
        help="Model enum names (multiple allowed with --payloads-only)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing eval outputs"
    )
    mode_group = parser.add_mutually_exclusive_group()
    for mode in (
        "archive",
        "archive-only",
        "validate-only",
        "payloads-only",
        "publish-parity",
    ):
        mode_group.add_argument(f"--{mode}", action="store_true")
    mode_group.add_argument(
        "--classify-yaml-changes",
        nargs="*",
        metavar="PATH",
        help="Print JSON payload-refresh classification for overlaid model YAMLs.",
    )
    parser.add_argument("--removed-yaml-paths", nargs="*", default=[])
    parser.add_argument("--full-roster", action="store_true")
    args = parser.parse_args(argv)

    if args.full_roster and not args.payloads_only:
        parser.error("--full-roster requires --payloads-only")
    if args.full_roster and args.models:
        parser.error("model-targeted refresh cannot use --full-roster")
    if args.publish_parity and args.models:
        parser.error("--publish-parity does not accept models")

    if args.classify_yaml_changes is not None:
        print(
            json.dumps(
                classify_yaml_changes(
                    args.classify_yaml_changes, args.removed_yaml_paths
                ),
                sort_keys=True,
            )
        )
        return 0

    checks = Checklist()
    models: list[Model] = args.models
    if len(models) > 1 and not args.payloads_only:
        parser.error("multiple models are supported only with --payloads-only")
    model = models[0] if models else None
    if args.payloads_only:
        if not models and not args.full_roster:
            parser.error("full payload refresh requires --full-roster")
        run_payload_refresh(checks, models)
    elif args.publish_parity:
        publish_parity_assets(checks)
    elif model is None:
        parser.error("model is required unless --payloads-only")
    elif (args.archive or args.archive_only) and not os.getenv("FIGSHARE_TOKEN"):
        parser.error("FIGSHARE_TOKEN must be set for archival")
    elif args.archive_only:
        run_archive(model, checks)
    else:
        banner(f"Ingesting model submission: {model.name}")
        energy_only = check_submission(
            model,
            checks,
            validate_runner=not args.validate_only,
        )
        overwrite = ("--overwrite",) if args.overwrite else ()
        run_model_steps(
            "STEP 2: Running evaluation scripts", EVAL_STEPS, model, checks,
            energy_only=energy_only, extra_args=("--auto-download", *overwrite),
        )  # fmt: skip
        run_model_steps(
            "STEP 3: Generating per-model figures", FIG_STEPS, model, checks,
            energy_only=energy_only,
        )  # fmt: skip
        if args.archive and not checks.n_failed:
            run_archive(model, checks)
            publish_parity_assets(checks)
        elif args.archive:
            checks.skip("Archival skipped because validation failed")
        if not args.validate_only and not checks.n_failed:
            run_payload_refresh(checks, (model,))

    banner("SUMMARY")
    print(checks.summary())
    if checks.n_failed:
        print(
            "\n⚠️  Some checks failed. Please review the PR checklist:\n"
            "   https://github.com/janosh/matbench-discovery/blob/main/"
            ".github/PULL_REQUEST_TEMPLATE/new_model.md"
        )
        return 1
    print("\n✅ All required checks passed!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
