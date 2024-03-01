from __future__ import annotations

from enum import StrEnum, unique
from typing import Self

from pymatviz.utils import styled_html_tag


class LabelEnum(StrEnum):
    """StrEnum with optional label and description attributes plus dict() method."""

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


@unique
class Key(LabelEnum):
    """Keys used to access dataframes columns."""

    arity = "arity", "Arity"
    bandgap_pbe = "bandgap_pbe", "PBE Band Gap"
    chem_sys = "chemical_system", "Chemical System"
    composition = "composition", "Composition"
    cse = "computed_structure_entry", "Computed Structure Entry"
    daf = "DAF", "Discovery Acceleration Factor"
    dft_energy = "uncorrected_energy", "DFT Energy"
    e_form = "e_form_per_atom_mp2020_corrected", "DFT E_form"
    e_form_pred = "e_form_per_atom_pred", "Predicted E_form"
    e_form_raw = "e_form_per_atom_uncorrected", "DFT E_form raw"
    e_form_wbm = "e_form_per_atom_wbm", "WBM E_form"
    each = "energy_above_hull", "E<sub>hull dist</sub>"
    each_pred = "e_above_hull_pred", "Predicted E<sub>hull dist</sub>"
    each_true = "e_above_hull_mp2020_corrected_ppd_mp", "E<sub>MP hull dist</sub>"
    each_wbm = "e_above_hull_wbm", "E<sub>WBM hull dist</sub>"
    final_struct = "relaxed_structure", "Relaxed Structure"
    forces = "forces", "Forces"
    form_energy = "formation_energy_per_atom", "Formation Energy (eV/atom)"
    formula = "formula", "Formula"
    init_struct = "initial_structure", "Initial Structure"
    magmoms = "magmoms", "Magnetic Moments"
    mat_id = "material_id", "Material ID"
    model_mean_each = "mean_pred_models", "Mean prediction all models"
    model_mean_err = "each_err_models", "Mean E<sub>hull dist</sub> error all models"
    model_std_each = "each_std_models", "Std. dev. all models"
    n_sites = "n_sites", "Number of Sites"
    site_nums = "site_nums", "Site Numbers", "Atomic numbers for each crystal site"
    spacegroup = "spacegroup", "Spacegroup"
    stress = "stress", "Stress"
    stress_trace = "stress_trace", "Stress Trace"
    struct = "structure", "Structure"
    task_id = "task_id", "Task ID"
    task_type = "task_type", "Task Type"
    train_task = "train_task", "Training Task"
    test_task = "test_task", "Test Task"
    targets = "targets", "Targets"
    # lowest WBM structures for a given prototype that isn't already in MP
    uniq_proto = "unique_prototype", "Unique Prototype"
    volume = "volume", "Volume (Å³)"
    wyckoff = "wyckoff_spglib", "Aflow-Wyckoff Label"  # relaxed structure Aflow label
    init_wyckoff = (
        "wyckoff_spglib_initial_structure",
        "Aflow-Wyckoff Label Initial Structure",
    )
    # number of structures in a model's training set
    train_set = "train_set", "Training Set"
    model_params = "model_params", "Model Params"  # model's parameter count
    model_type = "model_type", "Model Type"  # number of parameters in the model
    openness = "openness", "Openness"  # openness of data and code for a model


@unique
class Task(LabelEnum):
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
class Targets(LabelEnum):
    """Thermodynamic stability prediction task types."""

    E = "E", "energy"
    EFS = "EFS", "energy forces stress"
    EFSM = "EFSM", "energy forces stress magmoms"


@unique
class ModelType(LabelEnum):
    """Model types."""

    GNN = "GNN", "graph neural network"
    UIP = "UIP-GNN", "universal interatomic potential"
    BO_GNN = "BO-GNN", "GNN in a Bayesian optimization loop"
    Fingerprint = "Fingerprint", "models with structural fingerprint features"  # ex. RF
    Transformer = "Transformer", "transformer-based models"  # ex. Wrenformer
    RF = "RF", "random forest"


@unique
class Open(LabelEnum):
    """Openness of data and code for a model."""

    OSOD = "OSOD", "open source, open data"
    CSOD = "CSOD", "closed source, open data"
    OSCD = "OSCD", "open source, closed data"
    CSCD = "CSCD", "closed source, closed data"


@unique
class TestSubset(LabelEnum):
    """Test set subsets."""

    uniq_protos = "uniq_protos", "Unique Structure Prototypes"
    ten_k_most_stable = "10k_most_stable", "10k Most Stable"
    full = "full", "Full Test Set"


eV_per_atom = styled_html_tag(  # noqa: N816
    "(eV/atom)", tag="span", style="font-size: 0.8em; font-weight: lighter;"
)


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


class Model(LabelEnum):
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
