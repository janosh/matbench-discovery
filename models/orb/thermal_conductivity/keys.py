from enum import StrEnum, unique
from typing import Self


@unique
class Key(StrEnum):
    """Keys used to access dataframes columns, organized by semantic groups."""

    # Identifiers
    mat_id = "material_id"

    # metrics
    srme = "srme"
    sre = "sre"
    srd = "srd"

    # Thermal
    heat_capacity = "heat_capacity"

    # Phonon
    mode_weights = "mode_weights"
    has_imag_ph_modes = "has_imag_ph_modes"
    q_points = "q_points"
    ph_freqs = "ph_freqs"

    # Crystal Symmetry Properties
    spg_num = "spg_num"
    init_spg_num = "init_spg_num"
    final_spg_num = "final_spg_num"


class LabelEnum(StrEnum):
    """StrEnum with optional label and description attributes plus dict() methods."""

    def __new__(
        cls, val: str, label: str | None = None, desc: str | None = None
    ) -> Self:
        """Create a new class."""
        member = str.__new__(cls, val)
        member._value_ = val
        member.__dict__ |= dict(label=label, desc=desc)
        return member

    @property
    def label(self) -> str:
        """Make label read-only."""
        return self.__dict__["label"]

    @property
    def description(self) -> str:
        """Make description read-only."""
        return self.__dict__["desc"]


@unique
class MbdKey(LabelEnum):
    """Keys used to access dataframes columns."""

    # Thermal conductivity related keys
    kappa_tot_rta = (
        "kappa_tot_rta",
        "Total thermal conductivity (relaxation time approximation)",
    )
    kappa_tot_avg = "kappa_tot_avg", "Average total thermal conductivity"
    kappa_p_rta = "kappa_p_rta", "Particle-like thermal conductivity (RTA)"
    kappa_c = "kappa_c", "Thermal conductivity correction"
    mode_kappa_tot_rta = (
        "mode_kappa_tot_rta",
        "Mode-resolved total thermal conductivity (RTA)",
    )
    mode_kappa_tot_avg = (
        "mode_kappa_tot_avg",
        "Mode-resolved average thermal conductivity",
    )
    true_kappa_tot_avg = "true_kappa_tot_avg", "True average total thermal conductivity"
