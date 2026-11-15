"""Reproduce pair-residual selection using the author's frozen data/provenance."""

# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "pet-xs-pair-bo[train] @ git+https://github.com/nkwork9999/pet-xs-pair-bo.git@12ed3dce5231b46b7a81c3003d3b9cbcf5bab586",
# ]
# ///
from pet_pair_bo.training import main

if __name__ == "__main__":
    raise SystemExit(main())
