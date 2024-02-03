from __future__ import annotations

from enum import StrEnum, unique
from typing import Self

from pymatviz.utils import styled_html_tag


class DictableStrEnum(StrEnum):
    """StrEnum with optional description attributes and dict() method."""

    def __new__(cls, val: str, desc: str | None = None) -> Self:
        """Create a new class."""
        member = str.__new__(cls, val)
        member._value_ = val
        member.__dict__["desc"] = desc
        return member

    @property
    def description(self) -> str:
        """Make description read-only."""
        return self.__dict__["desc"]

    @classmethod
    def dict(cls) -> dict[str, str]:
        """Return the enum as a dictionary."""
        return {key: str(val) for key, val in cls.__members__.items()}


@unique
class Key(DictableStrEnum):
    """Keys used to access dataframes columns."""

    arity = "arity"
    bandgap_pbe = "bandgap_pbe"
    chem_sys = "chemical_system"
    composition = "composition"
    cse = "computed_structure_entry"
    dft_energy = "uncorrected_energy"
    e_form = "e_form_per_atom_mp2020_corrected"
    e_form_pred = "e_form_per_atom_pred"
    e_form_raw = "e_form_per_atom_uncorrected"
    e_form_wbm = "e_form_per_atom_wbm"
    each = "energy_above_hull"  # as returned by MP API
    each_pred = "e_above_hull_pred"
    each_true = "e_above_hull_mp2020_corrected_ppd_mp"
    each_wbm = "e_above_hull_wbm"
    final_struct = "relaxed_structure"
    forces = "forces"
    form_energy = "formation_energy_per_atom"
    formula = "formula"
    init_struct = "initial_structure"
    magmoms = "magmoms"
    mat_id = "material_id"
    model_mean_each = "Mean prediction all models"
    model_mean_err = "Mean error all models"
    model_std_each = "Std. dev. over models"
    n_sites = "n_sites"
    site_nums = "site_nums"
    spacegroup = "spacegroup"
    stress = "stress"
    stress_trace = "stress_trace"
    struct = "structure"
    task_id = "task_id"
    # lowest WBM structures for a given prototype that isn't already in MP
    uniq_proto = "unique_prototype"
    volume = "volume"
    wyckoff = "wyckoff_spglib"  # relaxed structure Aflow label
    init_wyckoff = "wyckoff_spglib_initial_structure"  # initial structure Aflow label


@unique
class Task(DictableStrEnum):
    """Thermodynamic stability prediction task types."""

    IS2RE = "IS2RE", "initial structure to relaxed energy"
    RS2RE = "RS2RE", "relaxed structure to relaxed energy"
    S2EFSM = "S2EFSM", "structure to energy force stress magmom"
    S2EFS = "S2EFS", "structure to energy force stress"
    # S2RE is for models that learned a discrete version of PES like CGCNN+P
    S2RE = "S2RE", "structure to relaxed energy"
    RP2RE = "RP2RE", "relaxed prototype to relaxed energy"
    IP2RE = "IP2RE", "initial prototype to relaxed energy"
    IS2E = "IS2E", "initial structure to energy"
    IS2RE_SR = "IS2RE-SR", "initial structure to relaxed energy after ML relaxation"


@unique
class Targets(DictableStrEnum):
    """Thermodynamic stability prediction task types."""

    E = "E", "energy"
    EFS = "EFS", "energy forces stress"
    EFSM = "EFSM", "energy forces stress magmoms"


@unique
class ModelType(DictableStrEnum):
    """Model types."""

    GNN = "GNN", "graph neural network"
    UIP = "UIP-GNN", "universal interatomic potential"
    BO_GNN = "BO-GNN", "GNN in a Bayesian optimization loop"
    Fingerprint = "Fingerprint", "models with structural fingerprint features"  # ex. RF
    Transformer = "Transformer", "transformer-based models"  # ex. Wrenformer
    RF = "RF", "random forest"


@unique
class Open(DictableStrEnum):
    """Openness of data and code for a model."""

    OSOD = "OSOD", "open source, open data"
    CSOD = "CSOD", "closed source, open data"
    OSCD = "OSCD", "open source, closed data"
    CSCD = "CSCD", "closed source, closed data"


ev_per_atom = styled_html_tag(
    "(eV/atom)", tag="span", style="font-size: 0.8em; font-weight: lighter;"
)


class Quantity(DictableStrEnum):
    """Quantity labels for plots."""

    n_atoms = "Atom Count"
    n_elems = "Element Count"
    crystal_sys = "Crystal system"
    spg_num = "Space group"
    n_wyckoff = "Number of Wyckoff positions"
    n_sites = "Number of atoms"
    energy_per_atom = f"Energy {ev_per_atom}"
    e_form = f"DFT E<sub>form</sub> {ev_per_atom}"
    e_above_hull = f"E<sub>hull dist</sub> {ev_per_atom}"
    e_above_hull_mp2020_corrected_ppd_mp = f"DFT E<sub>hull dist</sub> {ev_per_atom}"
    e_above_hull_pred = f"Predicted E<sub>hull dist</sub> {ev_per_atom}"
    e_above_hull_mp = f"E<sub>above MP hull</sub> {ev_per_atom}"
    e_above_hull_error = f"Error in E<sub>hull dist</sub> {ev_per_atom}"
    vol_diff = "Volume difference (A^3)"
    e_form_per_atom_mp2020_corrected = f"DFT E<sub>form</sub> {ev_per_atom}"
    e_form_per_atom_pred = f"Predicted E<sub>form</sub> {ev_per_atom}"
    material_id = "Material ID"
    band_gap = "Band gap (eV)"
    formula = "Formula"
    stress = "σ (eV/Å³)"  # noqa: RUF001
    stress_trace = "1/3 Tr(σ) (eV/Å³)"  # noqa: RUF001


class Model(DictableStrEnum):
    """Model labels for plots."""

    alignn_ff = "ALIGNN FF"
    alignn_pretrained = "ALIGNN Pretrained"
    alignn = "ALIGNN"
    bowsr_megnet = "BOWSR"
    cgcnn_p = "CGCNN+P"
    cgcnn = "CGCNN"
    chgnet_megnet = "CHGNet→MEGNet"
    chgnet = "CHGNet"
    dft = "DFT"
    gnome = "GNoME"
    m3gnet_direct = "M3GNet DIRECT"
    m3gnet_megnet = "M3GNet→MEGNet"
    m3gnet_ms = "M3GNet MS"
    m3gnet = "M3GNet"
    mace = "MACE"
    megnet_rs2re = "MEGNet RS2RE"
    megnet = "MEGNet"
    pfp = "PFP"
    voronoi_rf = "Voronoi RF"
    wbm = "WBM"
    wrenformer = "Wrenformer"
