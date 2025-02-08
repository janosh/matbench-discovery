"""Central argument parser for Matbench Discovery scripts."""

import argparse

from pymatviz.enums import Key

from matbench_discovery.data import Model
from matbench_discovery.enums import TestSubset

cli_parser = argparse.ArgumentParser(
    description="CLI flags for plotting and analysis scripts."
)
plot_group = cli_parser.add_argument_group(
    "plot", "Arguments for controlling figure generation"
)
plot_group.add_argument(
    "--models",
    nargs="*",
    type=Model,  # type: ignore[arg-type]
    choices=Model,
    default=list(Model),
    help="Models to analyze. If none specified, analyzes all models.",
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
