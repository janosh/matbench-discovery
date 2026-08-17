"""Central argument parser for Matbench Discovery scripts."""

import multiprocessing as mp
import os
import sys
from argparse import ArgumentParser, ArgumentTypeError

from matbench_discovery.enums import Model, TestSubset
from matbench_discovery.figs import PayloadMode


def parse_model(value: str) -> Model:
    """Parse a CLI model name into a Model enum member."""
    try:
        return Model.from_ref(value)
    except ValueError as exc:
        raise ArgumentTypeError(f"invalid model: {value}") from exc


def parse_key_migration(value: str) -> tuple[str, str]:
    """Parse one explicit ``OLD=NEW`` model-key migration."""
    old_key, separator, new_key = value.partition("=")
    if separator != "=" or not old_key or not new_key:
        raise ArgumentTypeError("model key migration must have the form OLD=NEW")
    return old_key, new_key


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
payload_mode_group = payload_group.add_mutually_exclusive_group()
payload_mode_group.add_argument(
    "--full-roster",
    action="store_true",
    help="Regenerate the full model roster without changing shared provenance.",
)
payload_mode_group.add_argument(
    "--migrate-provenance",
    action="store_true",
    help="Explicitly replace payload-wide computation provenance.",
)
payload_mode_group.add_argument(
    "--migrate-model-key",
    type=parse_key_migration,
    metavar="OLD=NEW",
    help="Explicitly migrate one immutable model key and its payload records.",
)
cli_args, _ignore_unknown = cli_parser.parse_known_args()
models_were_explicit = any(arg.partition("=")[0] == "--models" for arg in sys.argv)


def payload_mode() -> PayloadMode:
    """Return the one explicitly selected payload generation mode."""
    for selected, mode in (
        (cli_args.full_roster, PayloadMode.full_roster),
        (cli_args.migrate_provenance, PayloadMode.migrate_provenance),
        (cli_args.migrate_model_key, PayloadMode.migrate_model_key),
    ):
        if selected:
            if models_were_explicit:
                cli_parser.error(f"--{mode.value} cannot be combined with --models")
            return mode
    if models_were_explicit:
        return PayloadMode.targeted
    return cli_parser.error(
        "payload generation requires --models, --full-roster, "
        "--migrate-provenance, or --migrate-model-key"
    )


def is_full_model_run() -> bool:
    """Return whether the explicit payload mode covers the complete model roster."""
    return payload_mode() != PayloadMode.targeted


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
