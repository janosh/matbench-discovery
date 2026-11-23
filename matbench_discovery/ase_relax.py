"""Shared ASE optimizer and cell-filter resolution."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import FIRE, GoodOldQuasiNewton

if TYPE_CHECKING:
    from ase.filters import Filter
    from ase.optimize.optimize import Optimizer


# only the optimizers and cell filters benchmark submissions were validated with,
# keyed by their exact YAML names
OPTIMIZERS: dict[str, type[Optimizer]] = {"FIRE": FIRE, "GOQN": GoodOldQuasiNewton}
CELL_FILTERS: dict[str, type[Filter]] = {
    "FrechetCellFilter": FrechetCellFilter,
    "ExpCellFilter": ExpCellFilter,
}
# cell-filter name for a fixed cell, so CLI args (strings) can request one
NO_CELL_FILTER = "none"


def resolve_optimizer(name: str) -> type[Optimizer]:
    """Resolve and validate an ASE optimizer name."""
    try:
        return OPTIMIZERS[name]
    except KeyError:
        raise ValueError(
            f"Unknown ASE optimizer {name!r}, choose from {sorted(OPTIMIZERS)}"
        ) from None


def resolve_cell_filter(name: str | None) -> type[Filter] | None:
    """Resolve and validate an ASE cell-filter name, None or 'none' = fixed cell."""
    if name is None or name == NO_CELL_FILTER:
        return None
    try:
        return CELL_FILTERS[name]
    except KeyError:
        raise ValueError(
            f"Unknown ASE cell filter {name!r}, choose from "
            f"{[*sorted(CELL_FILTERS), NO_CELL_FILTER]}"
        ) from None


def canonical_optimizer_name(name: str) -> str:
    """Validate an optimizer name and return it for YAML/settings storage."""
    resolve_optimizer(name)
    return name


def canonical_filter_name(name: str | None) -> str | None:
    """Validate a cell-filter name, mapping 'none' to None (fixed cell)."""
    return None if resolve_cell_filter(name) is None else name
