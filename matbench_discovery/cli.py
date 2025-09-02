"""Central argument parser for Matbench Discovery scripts."""

import multiprocessing as mp
from argparse import ArgumentParser

from pymatviz.enums import Key

from matbench_discovery.enums import Model, TestSubset

cli_parser = ArgumentParser(
    description="CLI flags for eval, plot and analysis scripts."
)

cli_parser.add_argument(
    "--models",
    nargs="*",
    type=Model,  # type: ignore[arg-type]
    choices=Model,
    default=list(Model),
    help="Models to analyze. If none specified, analyzes all models.",
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
cli_parser.add_argument(
    "--timeout",
    type=int,
    default=30,
    help="Timeout in seconds for HTTP requests (default: 30).",
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
plot_group.add_argument(
    "--energy-type",
    type=str,
    default=Key.each,
    choices=[Key.e_form, Key.each],
    help="Whether to use formation energy or convex hull distance.",
)
plot_group.add_argument(
    "--show-non-compliant",
    action="store_true",
    help="Whether to show non-compliant models.",
)
plot_group.add_argument(
    "--use-full-rows",
    action="store_true",
    help="Whether to drop models that don't fit in complete rows.",
)
plot_group.add_argument(
    "--update-existing",
    action="store_true",
    help="Whether to update figures whose file paths already exist.",
)
cli_args, _ignore_unknown = cli_parser.parse_known_args()
