"""Enums used as keys/accessors for dicts and dataframes across Matbench Discovery."""

from enum import StrEnum, unique
from typing import Self

import pymatviz as pmv


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

    @classmethod
    def key_val_dict(cls) -> dict[str, str]:
        """Map of keys to values."""
        return {key: str(val) for key, val in cls.__members__.items()}

    @classmethod
    def val_label_dict(cls) -> dict[str, str | None]:
        """Map of values to labels."""
        return {str(val): val.label for val in cls.__members__.values()}

    @classmethod
    def val_desc_dict(cls) -> dict[str, str | None]:
        """Map of values to descriptions."""
        return {str(val): val.description for val in cls.__members__.values()}

    @classmethod
    def label_desc_dict(cls) -> dict[str | None, str | None]:
        """Map of labels to descriptions."""
        return {str(val.label): val.description for val in cls.__members__.values()}


eV_per_atom = pmv.html_tag(  # noqa: N816
    "(eV/atom)", tag="span", style="font-size: 0.8em; font-weight: lighter;"
)


@unique
class MbdKey(LabelEnum):
    """Keys used to access dataframes columns."""

    e_form_dft = (
        "e_form_per_atom_mp2020_corrected",
        f"DFT E<sub>form</sub> {eV_per_atom}",
    )
    e_form_raw = (
        "e_form_per_atom_uncorrected",
        f"DFT raw E<sub>form</sub> {eV_per_atom}",
    )
    e_form_wbm = "e_form_per_atom_wbm", f"WBM E<sub>form</sub> {eV_per_atom}"
    each_true = "e_above_hull_mp2020_corrected_ppd_mp", "E<sub>MP hull dist</sub>"
    each_mean_models = "each_mean_models", "E<sub>hull dist</sub> mean of models"
    each_err_models = "each_err_models", "E<sub>hull dist</sub> mean error of models"
    model_std_each = "each_std_models", "Std. dev. over models"
    init_wyckoff = (
        "wyckoff_spglib_initial_structure",
        "Aflow-Wyckoff Label Initial Structure",
    )
    openness = "openness", "Openness"
    dft_energy = "uncorrected_energy", "DFT Energy"
    each_wbm = "e_above_hull_wbm", "E<sub>WBM hull dist</sub>"


@unique
class Task(LabelEnum):
    """Thermodynamic stability prediction task types."""

    RP2RE = "RP2RE", "relaxed prototype to relaxed energy"
    RS2RE = "RS2RE", "relaxed structure to relaxed energy"
    S2E = "S2E", "structure to energy"
    # S2RE is for models that learned a discrete version of PES like CGCNN+P
    S2RE = "S2RE", "structure to relaxed energy"
    S2EF = "S2EF", "structure to energy, force"
    S2EFS = "S2EFS", "structure to energy, force, stress"
    S2EFSM = "S2EFSM", "structure to energy, force, stress, magmoms"
    IP2E = "IP2E", "initial prototype to energy"
    IS2E = "IS2E", "initial structure to energy"
    # IS2RE is for models that learned a discrete version of PES like CGCNN+P
    IS2RE = "IS2RE", "initial structure to relaxed energy"
    IS2RE_SR = (
        "IS2RE-SR",
        "initial structure to relaxed energy with structure relaxation",
    )
    IS2RE_BO = (
        "IS2RE-BO",
        "initial structure to relaxed energy with Bayesian optimization",
    )


@unique
class Targets(LabelEnum):
    """Thermodynamic stability prediction task types."""

    E = "E", "energy"
    EF = "EF", "energy forces"
    EFS = "EFS", "energy forces stress"
    EFSM = "EFSM", "energy forces stress magmoms"


@unique
class ModelType(LabelEnum):
    """Model types."""

    GNN = "GNN", "Graph Neural Network"
    UIP = "UIP", "Universal Interatomic Potential"
    BO_GNN = "BO-GNN", "GNN in a Bayesian Optimization Loop"
    Fingerprint = "Fingerprint", "Models with Structural Fingerprint Features"  # ex. RF
    Transformer = "Transformer", "Transformer-Based Models"  # ex. Wrenformer
    RF = "RF", "Random Forest"


@unique
class Open(LabelEnum):
    """Openness of data and code for a model."""

    OSOD = "OSOD", "open source, open data"
    CSOD = "CSOD", "closed source, open data"
    OSCD = "OSCD", "open source, closed data"
    CSCD = "CSCD", "closed source, closed data"


@unique
class TestSubset(LabelEnum):
    """Which subset of the test data to use for evaluation."""

    uniq_protos = "uniq_protos", "Unique Structure Prototypes"
    ten_k_most_stable = "10k_most_stable", "10k Most Stable"
    full = "full", "Full Test Set"


class Quantity(LabelEnum):
    """Quantity labels for plots."""

    n_atoms = "Atom Count"
    n_elems = "Element Count"
    crystal_sys = "Crystal system"
    spg_num = "Space group"
    n_wyckoff = "Number of Wyckoff positions"
    n_sites = "Number of atoms"
    energy_per_atom = f"Energy {eV_per_atom}"
    e_form = f"DFT E<sub>form</sub> {eV_per_atom}"
    e_above_hull = f"E<sub>hull dist</sub> {eV_per_atom}"
    e_above_hull_mp2020_corrected_ppd_mp = f"DFT E<sub>hull dist</sub> {eV_per_atom}"
    e_above_hull_pred = f"Predicted E<sub>hull dist</sub> {eV_per_atom}"
    e_above_hull_mp = f"E<sub>above MP hull</sub> {eV_per_atom}"
    e_above_hull_error = f"Error in E<sub>hull dist</sub> {eV_per_atom}"
    vol_diff = "Volume difference (A^3)"
    e_form_per_atom_mp2020_corrected = f"DFT E<sub>form</sub> {eV_per_atom}"
    e_form_per_atom_pred = f"Predicted E<sub>form</sub> {eV_per_atom}"
    material_id = "Material ID"
    band_gap = "Band gap (eV)"
    formula = "Formula"
    stress = "σ (eV/Å³)"  # noqa: RUF001
    stress_trace = "1/3 Tr(σ) (eV/Å³)"  # noqa: RUF001
