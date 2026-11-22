"""Geometry-relaxation settings: optimizer and cell-filter names resolved to matscikit.

Names follow ASE's class names so model YAMLs stay readable and stable; matscikit
implements the same algorithms (FIRE, GOQN; Frechet and Exp cell filters) with
geometry-identical trajectories.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from matscikit import RelaxResult, Structure
    from matscikit.model import Model as ForceModel
    from matscikit.options import BatchLimits, CellFilter, OnModelError, RelaxOptimizer

# canonical YAML name -> matscikit.relax(optimizer=...)
OPTIMIZERS: dict[str, RelaxOptimizer] = {"FIRE": "FIRE", "GOQN": "GOQN"}
OPTIMIZER_ALIASES = {"GoodOldQuasiNewton": "GOQN"}
# canonical YAML name -> matscikit.relax(cell_filter=...)
CELL_FILTERS: dict[str, CellFilter] = {
    "FrechetCellFilter": "frechet",
    "ExpCellFilter": "exp",
}
CELL_FILTER_ALIASES = {"frechet": "FrechetCellFilter", "exp": "ExpCellFilter"}
OPTIMIZER_LOOKUP = {key.casefold(): key for key in OPTIMIZERS} | {
    alias.casefold(): name for alias, name in OPTIMIZER_ALIASES.items()
}
CELL_FILTER_LOOKUP = {key.casefold(): key for key in CELL_FILTERS} | {
    alias.casefold(): name for alias, name in CELL_FILTER_ALIASES.items()
}


def canonical_optimizer_name(name: str) -> str:
    """Return the stable YAML name for an optimizer (case-insensitive, with aliases)."""
    try:
        return OPTIMIZER_LOOKUP[name.casefold()]
    except KeyError:
        raise ValueError(
            f"Unknown optimizer {name!r}, choose from {sorted(OPTIMIZERS)}"
        ) from None


def canonical_filter_name(name: str | None) -> str | None:
    """Return the stable YAML name for a cell filter, None (or 'none') = fixed cell."""
    if name is None or name.casefold() == "none":
        return None
    try:
        return CELL_FILTER_LOOKUP[name.casefold()]
    except KeyError:
        raise ValueError(
            f"Unknown cell filter {name!r}, choose from {sorted(CELL_FILTERS)}"
        ) from None


def relax_structures(
    structures: Sequence[Structure],
    model: ForceModel,
    *,
    optimizer: str,
    cell_filter: str | None,
    max_force: float,
    max_steps: int,
    fix_symmetry: float | None = None,
    cell_mask: Sequence[bool] | None = None,
    stress: bool = False,
    batch: BatchLimits | None = None,
    on_model_error: OnModelError = "raise",
    on_result: Callable[[int, RelaxResult], None] | None = None,
) -> list[RelaxResult]:
    """Relax structures in lock step (one batched model call per step).

    ``batch``, ``on_model_error`` and ``on_result`` pass through to ``matscikit.relax``
    (per-call structure/atom caps, isolating failing structures, streaming results).
    ``cell_mask`` (Voigt xx, yy, zz, yz, xz, xy) freezes masked cell strains, and
    ``stress=True`` reports the final stress of fixed-cell relaxations too. With
    ``fix_symmetry`` the start is idealized and forces/stress are symmetrized, like
    ASE's FixSymmetry.

    Returns one ``matscikit.RelaxResult`` per structure: relaxed ``structure``,
    ``converged``, ``n_steps`` and ``error`` (None, or why that structure failed, e.g.
    non-finite model output: its ``energy``/``forces``/``stress`` are then None), plus
    ``stress`` when the cell moves or ``stress=True``.
    """
    import matscikit

    filter_name = canonical_filter_name(cell_filter)
    return matscikit.relax(
        list(structures),
        model,
        f_max=max_force,
        max_steps=max_steps,
        optimizer=OPTIMIZERS[canonical_optimizer_name(optimizer)],
        cell_filter=CELL_FILTERS[filter_name] if filter_name else None,
        cell_mask=None if cell_mask is None or filter_name is None else list(cell_mask),
        fix_symmetry=fix_symmetry,
        stress=stress,
        batch=batch,
        on_model_error=on_model_error,
        on_result=on_result,
    )
