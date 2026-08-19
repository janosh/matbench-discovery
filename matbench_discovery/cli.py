"""Central argument parser for Matbench Discovery scripts."""

import multiprocessing as mp
import os
from argparse import ArgumentParser, ArgumentTypeError

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
)

cli_parser.add_argument(
    "--auto-download",
    action="store_true",
    help="Auto-confirm file downloads without prompting.",
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
    "--debug",
    type=int,
    default=0,
    help="If > 0, only analyze this many structures/items.",
)
cli_parser.add_argument(
    "--workers",
    type=int,
    default=max(1, mp.cpu_count() - 1),
    help="Number of processes to use for parallel tasks.",
)
cli_parser.add_argument(
    "--overwrite",
    action="store_true",
    help="Overwrite existing output files.",
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
cli_args, _ignore_unknown = cli_parser.parse_known_args()
models_were_explicit = cli_args.models is not models_arg.default


def is_full_model_run() -> bool:
    """Return whether payload generation explicitly selected the complete roster."""
    if cli_args.full_roster and models_were_explicit:
        cli_parser.error("--full-roster cannot be combined with --models")
    if not cli_args.full_roster and not models_were_explicit:
        cli_parser.error("payload generation requires --models or --full-roster")
    return cli_args.full_roster


def complete_models() -> list[Model]:
    """Return CLI-selected active models that have discovery prediction files."""
    return [
        model
        for model in cli_args.models
        if model.is_active and model.metrics.get("discovery", {}).get("pred_file")
    ]


# Set env var to auto-confirm file downloads when --auto-download is passed
if cli_args.auto_download:
    os.environ["MBD_AUTO_DOWNLOAD_FILES"] = "true"
