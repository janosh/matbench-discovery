"""Shared ASE optimizer and cell-filter resolution."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.filters import ExpCellFilter, FrechetCellFilter, UnitCellFilter
from ase.optimize import (
    BFGS,
    FIRE,
    FIRE2,
    LBFGS,
    BFGSLineSearch,
    GoodOldQuasiNewton,
    GPMin,
    LBFGSLineSearch,
)
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG

if TYPE_CHECKING:
    from ase.filters import Filter
    from ase.optimize.optimize import Optimizer


OPTIMIZERS: dict[str, type[Optimizer]] = {
    optimizer_class.__name__: optimizer_class
    for optimizer_class in (
        BFGS,
        BFGSLineSearch,
        FIRE,
        FIRE2,
        GPMin,
        GoodOldQuasiNewton,
        LBFGS,
        LBFGSLineSearch,
        SciPyFminBFGS,
        SciPyFminCG,
    )
}
OPTIMIZERS |= {"GOQN": GoodOldQuasiNewton, "QuasiNewton": BFGSLineSearch}
CELL_FILTERS: dict[str, type[Filter] | None] = {
    filter_class.__name__: filter_class
    for filter_class in (FrechetCellFilter, ExpCellFilter, UnitCellFilter)
}
CELL_FILTERS |= {"frechet": FrechetCellFilter, "exp": ExpCellFilter, "none": None}
OPTIMIZER_LOOKUP = {key.casefold(): value for key, value in OPTIMIZERS.items()}
CELL_FILTER_LOOKUP = {key.casefold(): value for key, value in CELL_FILTERS.items()}


def resolve_optimizer(name: str) -> type[Optimizer]:
    """Resolve and validate an ASE optimizer name case-insensitively."""
    try:
        return OPTIMIZER_LOOKUP[name.casefold()]
    except KeyError:
        raise ValueError(
            f"Unknown ASE optimizer {name!r}, choose from {sorted(OPTIMIZERS)}"
        ) from None


def resolve_cell_filter(name: str | None) -> type[Filter] | None:
    """Resolve and validate an ASE cell-filter name case-insensitively."""
    if name is None:
        return None
    try:
        return CELL_FILTER_LOOKUP[name.casefold()]
    except KeyError:
        raise ValueError(
            f"Unknown ASE cell filter {name!r}, choose from {sorted(CELL_FILTERS)}"
        ) from None


def canonical_optimizer_name(name: str) -> str:
    """Return the stable YAML name for an optimizer alias."""
    optimizer_class = resolve_optimizer(name)
    return "GOQN" if optimizer_class is GoodOldQuasiNewton else optimizer_class.__name__


def canonical_filter_name(name: str | None) -> str | None:
    """Return the stable YAML name for a cell-filter alias."""
    filter_class = resolve_cell_filter(name)
    return filter_class.__name__ if filter_class is not None else None
