"""Central argument parser for Matbench Discovery scripts."""

import multiprocessing as mp
import os
from argparse import ArgumentParser, ArgumentTypeError

import plotly.io as pio

from matbench_discovery.enums import Model, TestSubset


def parse_model(value: str) -> Model:
    """Parse a CLI model name into a Model enum member."""
    model = Model._missing_(value)
    if model is None:
        raise ArgumentTypeError(f"invalid model: {value}")
    return model


cli_parser = ArgumentParser(
    description="CLI flags for eval, plot and analysis scripts."
)

cli_parser.add_argument(
    "--auto-download",
    action="store_true",
    help="Auto-confirm file downloads without prompting.",
)

cli_parser.add_argument(
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
plot_group = cli_parser.add_argument_group(
    "plot", "Arguments for controlling figure generation"
)
plot_group.add_argument(
    "--test-subset",
    type=TestSubset,
    default=TestSubset.uniq_protos,
    choices=list(TestSubset),
    help="Which subset of the WBM test set to use for evaluation. "
    "Default is to only use unique Aflow protostructures. "
    "Training sets like MPtrj, sAlex and Omat24 were filtered to remove protostructures"
    " overlap with WBM, resulting in a slightly more out-of-distribution test set.",
)
cli_args, _ignore_unknown = cli_parser.parse_known_args()


def is_full_model_run() -> bool:
    """True when --models wasn't narrowed, i.e. the run covers all active models.
    Multi-model site figure payloads (site/src/figs) are only (over)written wholesale
    on full runs; filtered runs merge their entries into the committed payloads
    instead (see figs.write_site_payload), so they never clobber other models' data.
    """
    return set(cli_args.models) >= set(Model.active())


def complete_models() -> list[Model]:
    """Return CLI-selected models with complete discovery metrics."""
    return [model for model in cli_args.models if model.is_complete]


def shared_payload_test_subset() -> TestSubset:
    """Return the selected subset if it has one shared cohort across models."""
    if cli_args.test_subset == TestSubset.most_stable_10k:
        raise ValueError(
            "most_stable_10k is model-specific and cannot be represented in a shared "
            "multi-model payload"
        )
    return cli_args.test_subset


# Set env var to auto-confirm file downloads when --auto-download is passed
if cli_args.auto_download:
    os.environ["MBD_AUTO_DOWNLOAD_FILES"] = "true"

# Figures may open in browser tabs, but never steal focus.
for renderer_name in pio.renderers:
    renderer = pio.renderers[renderer_name]
    if hasattr(renderer, "autoraise"):
        renderer.autoraise = False
