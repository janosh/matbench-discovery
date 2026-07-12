"""Central dispatch for standard and backend-specific kappa adapters."""

from __future__ import annotations

from matbench_discovery.phonons.adapters.equflash import EquFlashKappaAdapter
from matbench_discovery.phonons.adapters.pet import PetKappaAdapter
from matbench_discovery.phonons.adapters.standard import (
    FairchemKappaAdapter,
    RelaxationResult,
    StandardKappaAdapter,
)


def get_kappa_adapter(model_key: str) -> StandardKappaAdapter:
    """Instantiate the custom adapter for a model or the standard ASE adapter."""
    if model_key == "pet_oam_xl_1_0_0":
        return PetKappaAdapter()
    if model_key.startswith("equflash"):
        return EquFlashKappaAdapter()
    if model_key.startswith(("eqv2_", "esen_", "equiformer_v3_")):
        return FairchemKappaAdapter()
    return StandardKappaAdapter()


__all__ = [
    "RelaxationResult",
    "StandardKappaAdapter",
    "get_kappa_adapter",
]
