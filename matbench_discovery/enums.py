"""Enums used as keys/accessors for dicts and dataframes across Matbench Discovery."""

import functools
import os
import re
from enum import EnumType, StrEnum, _EnumDict, auto, unique
from typing import Any, Self, TypeVar, cast

import pymatviz as pmv
import yaml

from matbench_discovery import DEFAULT_CACHE_DIR, PKG_DIR, ROOT
from matbench_discovery.remote.fetch import maybe_auto_download_file

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
    each_err_models = "each_err_models", "E<sub>hull dist</sub> mean error of models"
    openness = "openness", "Openness"

    init_protostructure_spglib = (
        "protostructure_spglib_initial_structure",
        "Protostructure Label for Initial Structure using spglib",
    )
    protostructure_spglib = (
        "protostructure_spglib",
        "Protostructure Label for Relaxed Structure using spglib",
    )
    international_spg_name = "international_spg_name", "International space group name"
    spg_num_diff = "spg_num_diff", "Difference in space group number"
    n_sym_ops_diff = "n_sym_ops_diff", "Difference in number of symmetry operations"
    structure_rmsd_vs_dft = "structure_rmsd_vs_dft", "RMSD of structure to DFT"

    # keep in sync with model-schema.yml
    missing_preds = "missing_preds", "Missing predictions"

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
    tortuosity = "tortuosity", "Tortuosity (unitless)"
    force_flips = "force_flips", "Force Flips (count)"
    energy_jump = "energy_jump", "Energy Jump (eV)"
    energy_diff_flips = "energy_diff_flips", "Energy Diff Flips (count)"
    force_total_variation = "force_total_variation", "Force Total Variation (eV/Å)"
    force_jump = "force_jump", "Force Jump (eV/Å)"
    pbe_wall_dist_mae = "pbe_wall_dist_mae", "PBE Wall Distance MAE (Å)"
    pbe_energy_mae = "pbe_energy_mae", "PBE Energy MAE (eV)"
    pbe_bond_length_error = "pbe_bond_length_error", "PBE Bond Length Error (Å)"
    pbe_well_depth_error = "pbe_well_depth_error", "PBE Well Depth Error (eV)"
    pbe_force_mae = "pbe_force_mae", "PBE Force MAE (eV/Å)"
    pbe_vib_freq_error = "pbe_vib_freq_error", "PBE Vib. Freq. Error (cm⁻¹)"


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
    EFSH_G = "EFSH_G", "EFSH<sub>G</sub>", "Energy with gradient-based Forces, Stress, and Hessian"  # fmt: skip


@unique
class ArchitectureType(LabelEnum):
    """Model architecture tags (GNN, transformer, classical ML, …)."""

    gnn = "gnn", "Graph Neural Network"
    transformer = "transformer", "Transformer"
    random_forest = "random_forest", "Random Forest"
    fingerprint = "fingerprint", "Structural Fingerprint"
    bayesian_optimization = "bayesian_optimization", "Bayesian Optimization"
    unknown = "unknown", "Unknown Architecture"


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


class MetaFiles(EnumType):
    """Metaclass of Files enum that adds base_dir class kwarg."""

    _base_dir: str

    def __new__(
        cls,
        name: str,
        bases: tuple[type, ...],
        namespace: _EnumDict,
        base_dir: str = DEFAULT_CACHE_DIR,
    ) -> "MetaFiles":
        """Create new Files enum with given base directory."""
        obj = super().__new__(cls, name, bases, namespace)
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
        # pass val to str.__new__ so the member's underlying string content is its
        # value. Without it, every member is the empty string "", so str.__eq__ and
        # str.__hash__ (inherited, not overridden) make all members compare equal and
        # hash to 0, silently collapsing any set/dict/`in` use of members.
        obj = str.__new__(cls, val)
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
    def url(self) -> str:
        """URL associated with the file."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement 'url' property"
        )

    @property
    def label(self) -> str:
        """Label associated with the file."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement 'label' property"
        )

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


# ruff: noqa: E501 (ignore long lines in class Model)
class Model(Files, base_dir=f"{ROOT}/models"):
    """Generated registry of non-aborted model YAML files.

    YAML ``model_key`` values determine member names. YAML metadata also provides
    lifecycle, metrics, hyperparams, package versions, links, and submission details.
    Regenerate the marked block with ``python scripts/generate_model_enum.py``.
    """

    # BEGIN GENERATED MODEL MEMBERS
    alchembert = auto(), "alchembert/alchembert.yml"
    alignn = auto(), "alignn/alignn.yml"
    allegro_mp_l_0_1 = auto(), "allegro/allegro-mp-l-0.1.yml"
    allegro_oam_l_0_1 = auto(), "allegro/allegro-oam-l-0.1.yml"
    alphanet_v1_mptrj = auto(), "alphanet/alphanet-v1-mptrj.yml"
    alphanet_v1_oam = auto(), "alphanet/alphanet-v1-oam.yml"
    bowsr = auto(), "bowsr/bowsr.yml"
    cgcnn = auto(), "cgcnn/cgcnn.yml"
    cgcnn_p = auto(), "cgcnn/cgcnn-p.yml"
    chgnet_0_3_0 = auto(), "chgnet/chgnet-0.3.0.yml"
    dpa3_v1_mptrj = auto(), "deepmd/dpa3-v1-mptrj.yml"
    dpa3_v1_openlam = auto(), "deepmd/dpa3-v1-openlam.yml"
    dpa3_v2_mptrj = auto(), "deepmd/dpa3-v2-mptrj.yml"
    dpa3_v2_openlam = auto(), "deepmd/dpa3-v2-openlam.yml"
    dpa_3_1_3m_ft = auto(), "deepmd/dpa-3.1-3m-ft.yml"
    dpa_3_1_mptrj = auto(), "deepmd/dpa-3.1-mptrj.yml"
    dpa_4_0_1_pro_mptrj = auto(), "deepmd/dpa-4.0.1-pro-mptrj.yml"
    dpa_4_0_pro_mptrj = auto(), "deepmd/dpa-4.0-pro-mptrj.yml"
    eqnorm_mptrj = auto(), "eqnorm/eqnorm-mptrj.yml"
    equflash_29m_oam = auto(), "equflash/equflash-29m-oam.yml"
    equflashv2_45m_oam = auto(), "equflash/equflashv2-45m-oam.yml"
    equiformer_v3_mp = auto(), "equiformer_v3/equiformer-v3-mp.yml"
    equiformer_v3_oam = auto(), "equiformer_v3/equiformer-v3-oam.yml"
    eqv2_m_omat_salex_mp = auto(), "eqv2/eqv2-m-omat-salex-mp.yml"
    eqv2_s_dens_mp = auto(), "eqv2/eqv2-s-dens-mp.yml"
    esen_30m_mp = auto(), "esen/esen-30m-mp.yml"
    esen_30m_oam = auto(), "esen/esen-30m-oam.yml"
    esnet = auto(), "esnet/esnet.yml"
    gnome = auto(), "gnome/gnome.yml"
    grace_1l_oam = auto(), "grace/grace-1l-oam.yml"
    grace_2l_mptrj = auto(), "grace/grace-2l-mptrj.yml"
    grace_2l_oam = auto(), "grace/grace-2l-oam.yml"
    grace_2l_oam_l = auto(), "grace/grace-2l-oam-l.yml"
    grace_3l_oam_l = auto(), "grace/grace-3l-oam-l.yml"
    hienet = auto(), "hienet/hienet.yml"
    m3gnet = auto(), "m3gnet/m3gnet.yml"
    mace_mp_0 = auto(), "mace/mace-mp-0.yml"
    mace_mpa_0 = auto(), "mace/mace-mpa-0.yml"
    matris_10m_mp = auto(), "matris/matris-10m-mp.yml"
    matris_10m_oam = auto(), "matris/matris-10m-oam.yml"
    matris_v050_mptrj = auto(), "matris/matris-v050-mptrj.yml"
    mattersim_v1_5m = auto(), "mattersim/mattersim-v1-5m.yml"
    megnet = auto(), "megnet/megnet.yml"
    nequip_mp_l_0_1 = auto(), "nequip/nequip-mp-l-0.1.yml"
    nequip_oam_l_0_1 = auto(), "nequip/nequip-oam-l-0.1.yml"
    nequip_oam_xl_0_1 = auto(), "nequip/nequip-oam-xl-0.1.yml"
    nequix_mp_1 = auto(), "nequix/nequix-mp-1.yml"
    nequix_mp_1_pft = auto(), "nequix/nequix-mp-1-pft.yml"
    orb_v2 = auto(), "orb/orb-v2.yml"
    orb_v2_mptrj = auto(), "orb/orb-v2-mptrj.yml"
    orb_v3 = auto(), "orb/orb-v3.yml"
    pet_oam_xl_1_0_0 = auto(), "pet/pet-oam-xl-1.0.0.yml"
    sevennet_0 = auto(), "sevennet/sevennet-0.yml"
    sevennet_l3i5 = auto(), "sevennet/sevennet-l3i5.yml"
    sevennet_mf_ompa = auto(), "sevennet/sevennet-mf-ompa.yml"
    sevennet_omni_i12 = auto(), "sevennet/sevennet-omni-i12.yml"
    tace_oam_l = auto(), "tace/tace-oam-l.yml"
    tace_oam_rra_preview = auto(), "tace/tace-oam-rra-preview.yml"
    tace_v1_oam_m = auto(), "tace/tace-v1-oam-m.yml"
    tece_oam_rra_1_0 = auto(), "tace/tece-oam-rra-1.0.yml"
    voronoi_rf = auto(), "voronoi_rf/voronoi-rf.yml"
    wrenformer = auto(), "wrenformer/wrenformer.yml"
    # END GENERATED MODEL MEMBERS

    @functools.cached_property  # cache to avoid re-reading same file multiple times
    def metadata(self) -> dict[str, Any]:
        """Metadata associated with the model."""
        yaml_path = f"{type(self).base_dir}/{self.rel_path}"
        with open(yaml_path, encoding="utf-8") as file:
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
        try:
            return self.metadata["pr_url"]
        except KeyError as exc:
            exc.add_note(f"{self.rel_path!r} missing required field 'pr_url'")
            raise

    @property
    def key(self) -> str:
        """Key associated with the file URL."""
        return self.metadata["model_key"]

    @property
    def family(self) -> str:
        """Model family directory derived from the YAML parent path."""
        return os.path.basename(os.path.dirname(self.rel_path))

    @property
    def metrics(self) -> dict[str, Any]:
        """Metrics associated with the model."""
        return self.metadata.get("metrics", {})

    @property
    def yaml_path(self) -> str:
        """YAML file path associated with the model."""
        return f"{type(self).base_dir}/{self.rel_path}"

    def _metric_pred_path(
        self,
        metrics_key: str,
        *,
        nested_key: str | None = None,
        required: bool = False,
    ) -> str | None:
        """Resolve metrics pred_file to a local path, downloading when a URL is set."""
        from matbench_discovery.data import file_ref_name, file_ref_url

        section = self.metrics.get(metrics_key)
        if not isinstance(section, dict):
            section = {}
        if nested_key is not None:
            section = section.get(nested_key) or {}
        pred_file = section.get("pred_file")
        if not (rel_path := file_ref_name(pred_file)):
            if required:
                raise ValueError(
                    f"metrics.{metrics_key}.pred_file not found in {self.rel_path!r}"
                )
            return None
        abs_path = f"{ROOT}/{rel_path}"
        if file_url := file_ref_url(pred_file):
            maybe_auto_download_file(file_url, abs_path, label=self.label)
        return abs_path

    @property
    def discovery_path(self) -> str:
        """Prediction file path associated with the model."""
        return cast("str", self._metric_pred_path("discovery", required=True))

    @property
    def geo_opt_path(self) -> str | None:
        """Geo-opt prediction path, downloading when a URL is present."""
        return self._metric_pred_path("geo_opt")

    @property
    def kappa_103_path(self) -> str | None:
        """Phonon kappa_103 prediction path, downloading when a URL is present."""
        return self._metric_pred_path("phonons", nested_key="kappa_103")

    @property
    def md_path(self) -> str | None:
        """MD prediction path, downloading when a URL is present."""
        return self._metric_pred_path("md")

    @property
    def is_active(self) -> bool:
        """Return whether the model participates in current benchmark views."""
        return self.metadata.get("lifecycle") == "active"

    @classmethod
    def active(cls) -> tuple[Self, ...]:
        """Active models included by default in plots, metrics, and eval scripts."""
        return tuple(model for model in cls if model.is_active)

    @classmethod
    def from_ref(cls, ref: str | Self) -> Self:
        """Resolve an enum name, value (case/dash-insensitive), key, or label to a
        member.

        Raises:
            ValueError: For unresolvable refs (with close-match suggestions).
        """
        if isinstance(ref, cls):
            return ref
        if member := cls.__members__.get(ref) or cls._missing_(ref):
            return member
        if member := next(
            (model for model in cls if ref in (model.key, model.label)), None
        ):
            return member
        return cls.from_label(ref)

    @classmethod
    def _missing_(cls, value: object) -> Self | None:
        """Normalize casing and punctuation before matching enum values.
        If no match is found, return None.

        This allows CLI arguments like --models mace-mp-0 to be recognized as mace_mp_0.
        """
        if isinstance(value, str):
            converted_value = re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")
            return cls._value2member_map_.get(converted_value)  # ty: ignore[invalid-return-type]

        return None


class DataFiles(Files):
    """Enum of data files with associated file directories and URLs."""

    mp_computed_structure_entries = (
        auto(),
        "mp/2023-02-07-mp-computed-structure-entries.json.gz",
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
        "wbm/2022-10-19-wbm-computed-structure-entries.jsonl.gz",
    )
    wbm_relaxed_atoms = auto(), "wbm/2024-08-04-wbm-relaxed-atoms.extxyz.zip"
    wbm_initial_structures = auto(), "wbm/2022-10-19-wbm-init-structs.jsonl.gz"
    wbm_initial_atoms = auto(), "wbm/2024-08-04-wbm-initial-atoms.extxyz.zip"
    wbm_summary = auto(), "wbm/2023-12-13-wbm-summary.csv.gz"
    phonondb_pbe_103_structures = (
        auto(),
        "phonons/2024-11-09-phononDB-PBE-103-structures.extxyz",
    )
    phonondb_pbe_103_kappa_no_nac = (
        auto(),
        "phonons/2024-11-09-kappas-phononDB-PBE-noNAC.json.gz",
    )
    wbm_dft_geo_opt_symprec_1e_2 = (
        auto(),
        "data/wbm/dft-geo-opt-symprec=1e-2-moyo=0.3.1.csv.gz",
    )
    wbm_dft_geo_opt_symprec_1e_5 = (
        auto(),
        "data/wbm/dft-geo-opt-symprec=1e-5-moyo=0.3.1.csv.gz",
    )
    dynamat_v1_0_md_trajectories = (
        auto(),
        "md/2026-06-29-dynamat-v1.0-reference-trajectories.h5",
    )

    @functools.cached_property
    def yaml(self) -> dict[str, dict[str, str]]:
        """YAML data associated with the file."""
        yaml_path = f"{PKG_DIR}/data-files.yml"

        with open(yaml_path, encoding="utf-8") as file:
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
            expected_md5 = self.yaml[key].get("md5")
            maybe_auto_download_file(self.url, abs_path, label=key, md5=expected_md5)
            if not os.path.isfile(abs_path):
                raise FileNotFoundError(
                    f"Failed to download and verify {key!r} from {self.url} "
                    f"to {abs_path} with {expected_md5=}."
                )
        return abs_path
