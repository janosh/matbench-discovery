"""Central argument parser for Matbench Discovery scripts."""

import sys
from argparse import ArgumentParser, ArgumentTypeError, Namespace
from collections.abc import Sequence

from matbench_discovery.enums import Model, TestSubset


def parse_model(value: str) -> Model:
    """Parse a CLI model name into a Model enum member."""
    try:
        return Model.from_ref(value)
    except ValueError as exc:
        raise ArgumentTypeError(f"invalid model: {value}") from exc


cli_parser = ArgumentParser(
    description="CLI flags for evaluation, payload and analysis scripts.",
    allow_abbrev=False,
    # scripts add their own flag groups, possibly again when tests re-import them
    conflict_handler="resolve",
)

models_arg = cli_parser.add_argument(
    "--models",
    nargs="+",
    type=parse_model,
    choices=Model,
    default=list(Model.active()),
    help="Models to analyze. If none specified, analyzes active models.",
)
cli_parser.add_argument(
    "-n",
    "--dry-run",
    action="store_true",
    help="Print what would be done without actually doing it.",
)
payload_group = cli_parser.add_argument_group(
    "payload", "Arguments for controlling payload generation"
)
payload_group.add_argument(
    "--test-subset",
    type=TestSubset,
    default=TestSubset.uniq_protos,
    choices=list(TestSubset),
    help="Which subset of the WBM test set to use for evaluation. "
    "Default is to only use unique Aflow protostructures. "
    "Training sets like MPtrj, sAlex and Omat24 were filtered to remove protostructures"
    " overlap with WBM, resulting in a slightly more out-of-distribution test set.",
)
payload_group.add_argument(
    "--full-roster",
    action="store_true",
    help="Regenerate the complete payload roster instead of selected model records.",
)


def parse_cli_args(raw_args: Sequence[str] | None = None) -> Namespace:
    """Strictly parse cli_parser (shared flags plus any groups a script added).

    Unknown flags exit with an error so typos can't silently fall back to defaults.
    The one exception is Jupyter: kernels pass their own flags (e.g. --f=kernel.json)
    when running a script's ``# %%`` cells interactively. The result is cached as
    ``cli_args``.
    """
    args, unknown = cli_parser.parse_known_args(raw_args)
    if unknown and "ipykernel" not in sys.modules:
        cli_parser.error(f"unrecognized arguments: {' '.join(unknown)}")
    globals()["cli_args"] = args
    return args


def _cli_args() -> Namespace:
    """Return the parsed CLI args, parsing sys.argv on first use."""
    return globals().get("cli_args") or parse_cli_args()


def __getattr__(name: str) -> Namespace:
    """Parse sys.argv lazily on first ``cli_args`` access (PEP 562), so importing
    this module never parses argv: scripts can first add their own argument groups
    to cli_parser, and their --help lists them.
    """
    if name == "cli_args":
        return _cli_args()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def is_full_model_run() -> bool:
    """Return whether payload generation explicitly selected the complete roster."""
    cli_args = _cli_args()
    models_were_explicit = cli_args.models is not models_arg.default
    if cli_args.full_roster == models_were_explicit:
        cli_parser.error("select exactly one of --models or --full-roster")
    return cli_args.full_roster


def complete_models() -> list[Model]:
    """Return CLI-selected active models that have discovery prediction files."""
    return [
        model
        for model in _cli_args().models
        if model.is_active and (model.metrics.get("discovery") or {}).get("pred_file")
    ]
