"""Enums used as keys/accessors for dicts and dataframes across Matbench Discovery."""

import abc
import builtins
import functools
import os
import sys
from enum import EnumMeta, StrEnum, _EnumDict, auto, unique
from typing import Any, Self, TypeVar

import plotly.express as px
import pymatviz as pmv
import yaml

from matbench_discovery import DEFAULT_CACHE_DIR, PKG_DIR, ROOT
from matbench_discovery.remote.fetch import download_file, maybe_auto_download_file

eV_per_atom = pmv.enums.eV_per_atom  # noqa: N816
T = TypeVar("T", bound="Files")


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
    e_above_hull_error = (
        "e_above_hull_error",
        f"Error in E<sub>hull dist</sub> {eV_per_atom}",
    )

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

    aflow_prototype = "aflow_prototype", "Aflow prototype"
    canonical_proto = "canonical_proto", "Canonical prototype"
    uniq_proto = "unique_prototype", "Unique prototype"

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

    # Diatomic curve metrics
    norm_auc = "norm_auc", "Norm. AUC (unitless)"
    smoothness = "smoothness", "Smoothness (eV/Å²)"
    tortuosity = "tortuosity", "Tortuosity (unitless)"
    force_flips = "force_flips", "Force Flips (count)"
    conservation = "conservation", "Conservation (eV/Å)"
    energy_jump = "energy_jump", "Energy Jump (eV)"
    energy_diff_flips = "energy_diff_flips", "Energy Diff Flips (count)"
    energy_grad_norm_max = "energy_grad_norm_max", "Energy Grad Norm Max (eV/Å)"
    force_total_variation = "force_total_variation", "Force Total Variation (eV/Å)"
    force_jump = "force_jump", "Force Jump (eV/Å)"
    energy_mae = "energy_mae", "Energy MAE vs Reference (eV)"
    force_mae = "force_mae", "Force MAE (eV/Å)"
    force_conservation = "force_conservation", "Force Conservation (eV/Å)"


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
    EFS_GM = "EFS_GM", "EFS<sub>G</sub>M", "Energy with gradient-based Forces and Stress; plus Magmoms"  # fmt: skip
    EFS_DM = "EFS_DM", "EFS<sub>D</sub>M", "Energy with direct Forces and Stress; plus Magmoms"  # fmt: skip


@unique
class ModelType(LabelEnum):
    """Model types."""

    GNN = "GNN", "Graph Neural Network"
    UIP = "UIP", "Universal Interatomic Potential"
    BO_GNN = "BO-GNN", "GNN in a Bayesian Optimization Loop"
    Fingerprint = "Fingerprint", "Models with Structural Fingerprint Features"  # ex. RF
    Transformer = "Transformer", "Transformer-Based Models"  # Wrenformer, AlchemBERT
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

    __test__ = False  # stop pytest from thinking this is a unit test

    uniq_protos = "unique_prototypes", "Unique Structure Prototypes"
    most_stable_10k = "most_stable_10k", "10k Most Stable Materials"
    full_test_set = "full_test_set", "Full Test Set"


class MetaFiles(EnumMeta):
    """Metaclass of Files enum that adds base_dir and (member|label)_map class
    properties.
    """

    _base_dir: str

    def __new__(
        cls,
        name: str,
        bases: tuple[type, ...],
        namespace: _EnumDict,
        base_dir: str = DEFAULT_CACHE_DIR,
        **kwargs: Any,
    ) -> "MetaFiles":
        """Create new Files enum with given base directory."""
        obj = super().__new__(cls, name, bases, namespace, **kwargs)
        obj._base_dir = base_dir
        return obj

    @property
    def base_dir(cls) -> str:
        """Base directory of the file."""
        return cls._base_dir


class Files(StrEnum, metaclass=MetaFiles):
    """Enum of data files with associated file directories and URLs."""

    def __new__(cls, val: str, file_path: str) -> Self:
        """Create a new member of the FileUrls enum with a given URL where to load the
        file from and directory where to save it to.
        """
        obj = str.__new__(cls)
        obj._value_ = val
        obj.__dict__ |= dict(file_path=file_path)
        return obj

    def __repr__(self) -> str:
        """String representation of the file."""
        return f"{type(self).__name__}.{self.name}"

    def __str__(self) -> str:
        """String representation of the file."""
        return self.name

    @property
    @abc.abstractmethod
    def url(self) -> str:
        """URL associated with the file."""

    @property
    @abc.abstractmethod
    def label(self) -> str:
        """Label associated with the file."""

    @property
    def rel_path(self) -> str:
        """Path of the file relative to the repo's ROOT directory."""
        return self.__dict__["file_path"]

    @classmethod
    def from_label(cls, label: str) -> Self:
        """Get enum member from pretty label."""
        file = next((attr for attr in cls if attr.label == label), None)
        if file is None:
            import difflib

            similar_labels = difflib.get_close_matches(label, [k.label for k in cls])
            raise ValueError(
                f"{label=} not found in {cls.__name__}. Did you mean one of {similar_labels}?"
            )
        return file


# ruff: noqa: E501, ERA001 (ignore long lines in class Model)
class Model(Files, base_dir=f"{ROOT}/models"):
    """Data files provided by Matbench Discovery.
    See https://janosh.github.io/matbench-discovery/contribute for data descriptions.
    """

    alchembert = auto(), "alchembert/alchembert.yml"

    # AlphaNet: https://arxiv.org/abs/2501.07155
    alphanet_mptrj = auto(), "alphanet/alphanet-mptrj.yml"

    # alignn with global pooling: https://arxiv.org/abs/2106.01829
    alignn = auto(), "alignn/alignn.yml"

    # alignn-ff with local pooling: https://arxiv.org/abs/2209.05554
    # Commented out because the model could not be evaluated due to OOM errors
    # see models/alignn_ff/readme.md
    # alignn_ff = auto(), "alignn/alignn-ff.yml"

    # BOWSR optimizer coupled with original megnet
    bowsr_megnet = auto(), "bowsr/bowsr.yml"

    # default CHGNet model from publication with 400,438 params
    chgnet_030 = auto(), "chgnet/chgnet-0.3.0.yml"

    # CGCNN 10-member ensemble
    cgcnn = auto(), "cgcnn/cgcnn.yml"

    # CGCNN 10-member ensemble with 5-fold training set perturbations
    cgcnn_p = auto(), "cgcnn/cgcnn+p.yml"

    # DeepMD-DPA3 models
    dpa3_v1_mptrj = auto(), "deepmd/dpa3-v1-mptrj.yml"
    dpa3_v1_openlam = auto(), "deepmd/dpa3-v1-openlam.yml"

    # FAIR-Chem
    eqv2_s_dens = auto(), "eqV2/eqV2-s-dens-mp.yml"
    eqv2_m = auto(), "eqV2/eqV2-m-omat-salex-mp.yml"

    # GRACE: https://arxiv.org/abs/2311.16326v2
    grace_2l_mptrj = auto(), "grace/grace-2l-mptrj.yml"
    grace_2l_oam = auto(), "grace/grace-2l-oam.yml"
    grace_1l_oam = auto(), "grace/grace-1l-oam.yml"

    # GNoME - Nequip architecture trained on Google's proprietary data. Weights
    # are not publicly available and so these results cannot be reproduced.
    gnome = auto(), "gnome/gnome.yml"

    # original M3GNet straight from publication, not re-trained
    m3gnet_ms = auto(), "m3gnet/m3gnet.yml"
    # m3gnet_direct = auto(), "M3GNet DIRECT"
    # m3gnet_ms = auto(), "M3GNet MS"

    # MACE-MP-0 medium as published in https://arxiv.org/abs/2401.00096 trained on MPtrj
    mace_mp_0 = auto(), "mace/mace-mp-0.yml"
    mace_mpa_0 = auto(), "mace/mace-mpa-0.yml"  # trained on MPtrj and Alexandria

    # MatRIS-v0.5.0-MPtrj
    matris_v050_mptrj = auto(), "matris/matris-v050-mptrj.yml"
    
    # MatterSim - M3gNet architecture trained on propertary MSFT data. Weights
    # are open-sourced.
    mattersim_v1_5m = auto(), "mattersim/mattersim-v1-5m.yml"

    # original MEGNet straight from publication, not re-trained
    megnet = auto(), "megnet/megnet.yml"

    # ORB
    orb = auto(), "orb/orb.yml"
    orb_mptrj = auto(), "orb/orb-mptrj.yml"

    # SevenNet trained on MPtrj
    sevennet_0 = auto(), "sevennet/sevennet-0.yml"
    sevennet_l3i5 = auto(), "sevennet/sevennet-l3i5.yml"

    # Magpie composition+Voronoi tessellation structure features + sklearn random forest
    voronoi_rf = auto(), "voronoi_rf/voronoi-rf.yml"

    # wrenformer 10-member ensemble
    wrenformer = auto(), "wrenformer/wrenformer.yml"

    # --- Model Combos
    # # CHGNet-relaxed structures fed into MEGNet for formation energy prediction
    # chgnet_megnet = "chgnet/2023-03-06-chgnet-0.2.0-wbm-IS2RE.csv.gz"
    # # M3GNet-relaxed structures fed into MEGNet for formation energy prediction
    # m3gnet_megnet = "m3gnet/2022-10-31-m3gnet-wbm-IS2RE.csv.gz"
    # megnet_rs2re = "megnet/2023-08-23-megnet-wbm-RS2RE.csv.gz"

    @functools.cached_property  # cache to avoid re-reading same file multiple times
    def metadata(self) -> dict[str, Any]:
        """Metadata associated with the model."""
        yaml_path = f"{type(self).base_dir}/{self.rel_path}"
        with open(yaml_path) as file:
            data = yaml.safe_load(file)

        if not isinstance(data, dict):
            raise TypeError(f"{yaml_path!r} does not contain valid YAML metadata")

        return data

    @property
    def label(self) -> str:
        """Pretty label associated with the model."""
        return self.metadata["model_name"]

    @property
    def pr_url(self) -> str:
        """Pull request URL in which the model was originally added to the repo."""
        return self.metadata["pr_url"]

    @property
    def key(self) -> str:
        """Key associated with the file URL."""
        return self.metadata["model_key"]

    @property
    def metrics(self) -> dict[str, Any]:
        """Metrics associated with the model."""
        return self.metadata.get("metrics", {})

    @property
    def yaml_path(self) -> str:
        """YAML file path associated with the model."""
        return f"{type(self).base_dir}/{self.rel_path}"

    @property
    def discovery_path(self) -> str:
        """Prediction file path associated with the model."""
        rel_path = self.metrics.get("discovery", {}).get("pred_file")
        file_url = self.metrics.get("discovery", {}).get("pred_file_url")
        if not rel_path:
            raise ValueError(
                f"metrics.discovery.pred_file not found in {self.rel_path!r}"
            )
        abs_path = f"{ROOT}/{rel_path}"
        maybe_auto_download_file(file_url, abs_path, label=self.label)
        return abs_path

    @property
    def geo_opt_path(self) -> str | None:
        """File path associated with the file URL if it exists, otherwise
        download the file first, then return the path.
        """
        geo_opt_metrics = self.metrics.get("geo_opt", {})
        if geo_opt_metrics in ("not available", "not applicable"):
            return None
        rel_path = geo_opt_metrics.get("pred_file")
        file_url = geo_opt_metrics.get("pred_file_url")
        if not rel_path:
            raise ValueError(
                f"metrics.geo_opt.pred_file not found in {self.rel_path!r}"
            )
        abs_path = f"{ROOT}/{rel_path}"
        maybe_auto_download_file(file_url, abs_path, label=self.label)
        return abs_path

    @property
    def kappa_103_path(self) -> str | None:
        """File path associated with the file URL if it exists, otherwise
        download the file first, then return the path.
        """
        phonons_metrics = self.metrics.get("phonons", {})
        if phonons_metrics in ("not available", "not applicable"):
            return None
        rel_path = phonons_metrics.get("kappa_103", {}).get("pred_file")
        file_url = phonons_metrics.get("kappa_103", {}).get("pred_file_url")
        if not rel_path:
            raise ValueError(
                f"metrics.phonons.kappa_103.pred_file not found in {self.rel_path!r}"
            )
        abs_path = f"{ROOT}/{rel_path}"
        maybe_auto_download_file(file_url, abs_path, label=self.label)
        return abs_path


class DataFiles(Files):
    """Enum of data files with associated file directories and URLs."""

    mp_computed_structure_entries = (
        auto(),
        ("mp/2023-02-07-mp-computed-structure-entries.json.gz"),
    )
    mp_elemental_ref_entries = (
        auto(),
        "mp/2023-02-07-mp-elemental-reference-entries.json.gz",
    )
    # this file was originally generated on 2023-01-10, but was updated on 2025-02-01
    # to include moyopy-powered symmetry analysis of MP ground state structures
    mp_energies = auto(), "mp/2025-02-01-mp-energies.csv.gz"
    mp_patched_phase_diagram = auto(), "mp/2023-02-07-ppd-mp.pkl.gz"
    mp_trj_json_gz = auto(), "mp/2022-09-16-mp-trj.json.gz"
    mp_trj_extxyz = auto(), "mp/2024-09-03-mp-trj.extxyz.zip"
    # snapshot of every task (calculation) in MP as of 2023-03-16 (14 GB)
    all_mp_tasks = auto(), "mp/2023-03-16-all-mp-tasks.zip"

    wbm_computed_structure_entries = (
        auto(),
        ("wbm/2022-10-19-wbm-computed-structure-entries.json.bz2"),
    )
    wbm_relaxed_atoms = auto(), "wbm/2024-08-04-wbm-relaxed-atoms.extxyz.zip"
    wbm_initial_structures = auto(), "wbm/2022-10-19-wbm-init-structs.json.bz2"
    wbm_initial_atoms = auto(), "wbm/2024-08-04-wbm-initial-atoms.extxyz.zip"
    wbm_cses_plus_init_structs = (
        auto(),
        ("wbm/2022-10-19-wbm-computed-structure-entries+init-structs.json.bz2"),
    )
    wbm_summary = auto(), "wbm/2023-12-13-wbm-summary.csv.gz"
    alignn_checkpoint = auto(), "2023-06-02-pbenner-best-alignn-model.pth.zip"
    phonondb_pbe_103_structures = (
        auto(),
        ("phonons/2024-11-09-phononDB-PBE-103-structures.extxyz"),
    )
    phonondb_pbe_103_kappa_no_nac = (
        auto(),
        ("phonons/2024-11-09-kappas-phononDB-PBE-noNAC.json.gz"),
    )
    wbm_dft_geo_opt_symprec_1e_2 = (
        auto(),
        "data/wbm/dft-geo-opt-symprec=1e-2-moyo=0.3.1.csv.gz",
    )
    wbm_dft_geo_opt_symprec_1e_5 = (
        auto(),
        "data/wbm/dft-geo-opt-symprec=1e-5-moyo=0.3.1.csv.gz",
    )

    @functools.cached_property
    def yaml(self) -> dict[str, dict[str, str]]:
        """YAML data associated with the file."""
        yaml_path = f"{PKG_DIR}/data-files.yml"

        with open(yaml_path) as file:
            yaml_data = yaml.safe_load(file)

        if self.name not in yaml_data:
            raise ValueError(f"{self.name=} not found in {yaml_path}")

        return yaml_data

    @property
    def url(self) -> str:
        """URL associated with the file."""
        url = self.yaml[self.name].get("url")
        if url is None:
            raise ValueError(f"{self.name!r} does not have a URL")
        return url

    @property
    def label(self) -> str:
        """No pretty label for DataFiles, use name instead."""
        return self.name

    @property
    def description(self) -> str:
        """Description associated with the file."""
        return self.yaml[self.name]["description"]

    @property
    def path(self) -> str:
        """File path associated with the file URL if it exists, otherwise
        download the file first, then return the path.
        """
        key, rel_path = self.name, self.rel_path

        if rel_path not in self.yaml[key]["path"]:
            raise ValueError(f"{rel_path=} does not match {self.yaml[key]['path']}")

        abs_path = f"{type(self).base_dir}/{rel_path}"
        if not os.path.isfile(abs_path):
            is_ipython = hasattr(builtins, "__IPYTHON__")
            # default to 'y' if not in interactive session, and user can't answer
            answer = (
                input(
                    f"{abs_path!r} associated with {key=} does not exist. Download it "
                    "now? This will cache the file for future use. [y/n] "
                )
                if is_ipython or sys.stdin.isatty()
                else "y"
            )
            if answer.lower().strip() == "y":
                print(f"Downloading {key!r} from {self.url} to {abs_path}")
                download_file(abs_path, self.url)
        return abs_path


# register pretty labels to use instead of enum keys in plotly axes and legends
px.defaults.labels |= {key.name: key.label for key in (*Model, *MbdKey, *pmv.enums.Key)}
