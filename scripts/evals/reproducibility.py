"""Compare independent reruns of a task to its registered pred_file and write the
verdict to metrics.<task>.reproducibility in the model YAML (rules in
matbench_discovery/metrics/reproducibility.py, usage in contributing.md).

Example:
    uv run scripts/evals/reproducibility.py --task diatomics --models mace_mp_0 \
        --run "published run (2026-06-28, H100)" \
        models/mace/mace-mp-0/2026-06-28-diatomics.json.gz [URL]
"""

import os

from matbench_discovery import ROOT, file_digest, repo_relative_path
from matbench_discovery.cli import cli_parser, parse_cli_args
from matbench_discovery.data import (
    make_file_ref,
    update_yaml_file,
    validate_artifact_path,
)
from matbench_discovery.enums import Model, resolve_verified_file
from matbench_discovery.metrics import reproducibility as repro


def main() -> int:
    """Score the registered file and every run, then write the verdict block."""
    # flags added here, not at import, so --task never becomes required elsewhere
    group = cli_parser.add_argument_group("reproducibility")
    group.add_argument("--task", choices=("diatomics", "md"), required=True)
    group.add_argument(
        "--run",
        nargs="+",
        action="append",
        required=True,
        metavar="LABEL PATH [URL]",
        help="Independent run (repeatable), PATH downloaded from URL if missing",
    )
    group.add_argument(
        "--repeat",
        nargs="+",
        default=[],
        metavar="PATH",
        help=">= 2 identical same-env runs to measure the noise floor",
    )
    group.add_argument("--cause", choices=repro.REPRO_CAUSES)
    group.add_argument("--note", help="One sentence, required unless reproduced")
    args = parse_cli_args()
    if len(args.models) != 1:
        cli_parser.error("run paths belong to one model: pass exactly one --models")
    if any(len(run_spec) not in (2, 3) for run_spec in args.run):
        cli_parser.error(f"--run takes LABEL PATH [URL], got {args.run}")
    if len(args.repeat) == 1:
        cli_parser.error("--repeat takes at least 2 paths")
    model: Model = args.models[0]
    task: repro.ReproTask = args.task
    if (registered_path := getattr(model, f"{task}_path")) is None:
        raise ValueError(f"{model.name} has no registered metrics.{task}.pred_file")

    runs: list[dict[str, object]] = []
    run_paths: list[str] = []
    for label, path, *maybe_url in args.run:
        url = maybe_url[0] if maybe_url else None
        rel_path = repo_relative_path(path)
        validate_artifact_path(rel_path, f"models/{model.family}/{model.key}")
        abs_path = resolve_verified_file(
            f"{ROOT}/{rel_path}", url=url, label=label, md5=None
        )
        file_ref = make_file_ref(
            rel_path,
            url=url,
            size=os.path.getsize(abs_path),
            md5=file_digest(abs_path, "md5"),
        )
        runs.append({"label": label, "file": file_ref})
        run_paths.append(abs_path)

    repeat_paths = [f"{ROOT}/{repo_relative_path(path)}" for path in args.repeat]
    registered, *scores = repro.score_files(
        task, model, [registered_path, *run_paths, *repeat_paths]
    )
    block = repro.make_block(
        task,
        registered,
        scores[: len(runs)],
        runs=runs,
        repeat_scores=scores[len(runs) :],
        cause=args.cause,
        note=args.note,
    )
    print(f"{model.name} {task}: {block}")
    if not args.dry_run:
        key = f"metrics.{task}.reproducibility"
        update_yaml_file(model.yaml_path, key, block, preserve_existing=False)
        print(f"Wrote {key} to {model.yaml_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
