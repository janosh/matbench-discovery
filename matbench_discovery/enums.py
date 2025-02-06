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


eV_per_atom = pmv.html_tag(  # noqa: N816
    "(eV/atom)", tag="span", style="font-size: 0.8em; font-weight: lighter;"
)


@unique
class MbdKey(LabelEnum):
    """Keys used to access dataframes columns."""

    dft_energy = "uncorrected_energy", "DFT Energy"
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
    each_wbm = "e_above_hull_wbm", "E<sub>WBM hull dist</sub>"
    each_mean_models = "each_mean_models", "E<sub>hull dist</sub> mean of models"
    each_err_models = "each_err_models", "E<sub>hull dist</sub> mean error of models"
    model_std_each = "each_std_models", "Std. dev. over models"
    openness = "openness", "Openness"
    e_above_hull_error = f"Error in E<sub>hull dist</sub> {eV_per_atom}"

    init_wyckoff_spglib = (
        "wyckoff_spglib_initial_structure",
        "Protostructure Label for Initial Structure using spglib",
    )
    wyckoff_spglib = (
        "wyckoff_spglib",
        "Protostructure Label for Relaxed Structure using spglib",
    )
    international_spg_name = "international_spg_name", "International space group name"
    spg_num_diff = "spg_num_diff", "Difference in space group number"
    n_sym_ops_diff = "n_sym_ops_diff", "Difference in number of symmetry operations"
    structure_rmsd_vs_dft = "structure_rmsd_vs_dft", "RMSD of structure to DFT"
    sym_prop = "symmetry_property", "Symmetry property"

    # keep in sync with model-schema.yml
    missing_preds = "missing_preds", "Missing predictions"
    missing_percent = "missing_percent", "Missing predictions (percent)"

    aflow_prototype = "aflow_prototype"
    canonical_proto = "canonical_proto"
    uniq_proto = "unique_prototype"

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

    E = "E", "Energy"
    EF_G = "EF_G", "EF<sub>G</sub>", "Energy with gradient-based Forces"
    EF_D = "EF_D", "EF<sub>D</sub>", "Energy with direct Forces"
    EFS_G = "EFS_G", "EFS<sub>G</sub>", "Energy with gradient-based Forces and Stress"
    EFS_D = "EFS_D", "EFS<sub>D</sub>", "Energy with direct Forces and Stress"
    EFS_GM = "EFS_GM", "EFS<sub>G</sub>M", "Energy with gradient-based Forces and Stress; plus Magmoms"  # fmt: skip  # noqa: E501
    EFS_DM = "EFS_DM", "EFS<sub>D</sub>M", "Energy with direct Forces and Stress; plus Magmoms"  # fmt: skip  # noqa: E501


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

    uniq_protos = "unique_prototypes", "Unique Structure Prototypes"
    most_stable_10k = "most_stable_10k", "10k Most Stable Materials"
    full_test_set = "full_test_set", "Full Test Set"
