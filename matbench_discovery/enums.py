"""Enums used as keys/accessors for dicts and dataframes across Matbench Discovery."""

import functools
import keyword
import os
import re
from enum import EnumType, StrEnum, _EnumDict, auto, unique
from glob import glob
from typing import Any, Self, cast

import pymatviz as pmv
import yaml

from matbench_discovery import DEFAULT_CACHE_DIR, PKG_DIR, ROOT, file_digest
from matbench_discovery.remote.fetch import maybe_auto_download_file

eV_per_atom = pmv.enums.eV_per_atom  # noqa: N816


def resolve_verified_file(
    abs_path: str,
    *,
    url: str | None,
    label: str,
    md5: str | None,
) -> str:
    """Resolve a cached file and reject it unless its declared MD5 verifies."""
    # verified cache hit: return before maybe_auto_download_file hashes the file again
    if (
        md5 is not None
        and os.path.isfile(abs_path)
        and file_digest(abs_path, "md5") == md5
    ):
        return abs_path
    if url is not None:
        maybe_auto_download_file(url, abs_path, label=label, md5=md5)
    if not os.path.isfile(abs_path):
        raise FileNotFoundError(
            f"Failed to resolve {label!r} at {abs_path} with {url=}"
        )
    if md5 is not None:
        actual_md5 = file_digest(abs_path, "md5")
        if actual_md5 != md5:
            raise ValueError(
                f"Stale cached file for {label!r} at {abs_path}: "
                f"expected MD5 {md5}, got {actual_md5}"
            )
    return abs_path


class LabelEnum(StrEnum):
    """StrEnum whose members carry an optional pretty label."""

    def __new__(cls, val: str, label: str | None = None) -> Self:
        """Create a new class."""
        member = str.__new__(cls, val)
        member._value_ = val
        member.__dict__["label"] = label
        return member

    @property
    def label(self) -> str:
        """Make label read-only."""
        return self.__dict__["label"]


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
    EF_G = "EF_G", "EF<sub>G</sub>"  # Energy with gradient-based Forces
    EF_D = "EF_D", "EF<sub>D</sub>"  # Energy with direct Forces
    EFS_G = "EFS_G", "EFS<sub>G</sub>"  # gradient-based Forces and Stress
    EFS_D = "EFS_D", "EFS<sub>D</sub>"  # direct Forces and Stress
    EFS_GM = "EFS_GM", "EFS<sub>G</sub>M"  # gradient-based F and S plus Magmoms
    EFS_DM = "EFS_DM", "EFS<sub>D</sub>M"  # direct F and S plus Magmoms
    EFSH_G = "EFSH_G", "EFSH<sub>G</sub>"  # gradient-based F, S and Hessian


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
    ) -> MetaFiles:
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
        """Create a member identified by ``val`` and holding its relative file path."""
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


MODEL_KEY_PATTERN = re.compile(r"[a-z0-9]+(?:[.-][a-z0-9]+)*")


def _model_members() -> dict[str, tuple[auto, str]]:
    """Map each non-aborted model YAML's enum name (its model_key with - and . as _)
    to its path relative to models/, rejecting invalid or colliding keys.
    """
    members: dict[str, tuple[auto, str]] = {}
    key_owners: dict[str, str] = {}
    for yaml_path in sorted(glob(f"{ROOT}/models/[!_]*/[!_]*.yml")):
        with open(yaml_path, encoding="utf-8") as file:
            # libyaml parses all model YAMLs ~8x faster than PyYAML's pure-Python loader
            metadata = yaml.load(file, Loader=yaml.CSafeLoader)
        if not isinstance(metadata, dict):
            raise TypeError(f"{yaml_path} must contain a YAML mapping")
        model_key = metadata.get("model_key")
        if not isinstance(model_key, str) or not MODEL_KEY_PATTERN.fullmatch(model_key):
            raise TypeError(f"{yaml_path} has invalid {model_key=}")
        if owner := key_owners.get(model_key):
            raise ValueError(f"{model_key=} collides in {owner} and {yaml_path}")
        key_owners[model_key] = yaml_path
        if metadata.get("lifecycle") == "aborted":
            continue
        name = model_key.replace("-", "_").replace(".", "_")
        if not name.isidentifier() or keyword.iskeyword(name) or name in members:
            raise ValueError(f"{yaml_path} has invalid or duplicate enum name {name!r}")
        rel_path = os.path.relpath(yaml_path, f"{ROOT}/models").replace("\\", "/")
        members[name] = auto(), rel_path
    return dict(sorted(members.items()))  # name order, as the members were written


class Model(Files, base_dir=f"{ROOT}/models"):
    """Registry of non-aborted model YAML files, built from models/ at import.

    YAML ``model_key`` values determine member names. YAML metadata also provides
    lifecycle, metrics, hyperparams, package versions, links, and submission details.
    """

    locals().update(_model_members())  # one member per non-aborted model YAML

    @functools.cached_property  # cache to avoid re-reading same file multiple times
    def metadata(self) -> dict[str, Any]:
        """Metadata associated with the model."""
        with open(self.yaml_path, encoding="utf-8") as file:
            data = yaml.safe_load(file)

        if not isinstance(data, dict):
            raise TypeError(f"{self.yaml_path!r} does not contain valid YAML metadata")

        return data

    @property
    def label(self) -> str:
        """Pretty label associated with the model."""
        return self.metadata["model_name"]

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
        return self.metadata.get("metrics") or {}  # bare `metrics:` loads as None

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
        file_key: str = "pred_file",
    ) -> str | None:
        """Resolve a metrics file ref to a local path, downloading when a URL is set."""
        from matbench_discovery.data import (
            file_ref_name,
            file_ref_url,
            validate_artifact_path,
        )

        section = self.metrics.get(metrics_key) or {}
        for key in nested_key.split(".") if nested_key else ():
            section = section.get(key) or {}
        pred_file = section.get(file_key) or {}
        if not (rel_path := file_ref_name(pred_file)):
            if required:
                raise ValueError(
                    f"metrics.{metrics_key}.{file_key} not found in {self.rel_path!r}"
                )
            return None
        validate_artifact_path(rel_path, f"models/{self.family}/{self.key}")
        abs_path = f"{ROOT}/{rel_path}"
        return resolve_verified_file(
            abs_path,
            url=file_ref_url(pred_file),
            label=self.label,
            md5=pred_file.get("md5"),
        )

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
    def kappa_103_phonon_path(self) -> str | None:
        """Harmonic phonon sidecar (band path, eigenvectors, DOS) from
        ``kappa_103.harmonic`` for a harmonic-only rerun, else from ``kappa_103``.
        """
        kappa_103 = (self.metrics.get("phonons") or {}).get("kappa_103") or {}
        nested_key = "kappa_103.harmonic" if "harmonic" in kappa_103 else "kappa_103"
        return self._metric_pred_path(
            "phonons", nested_key=nested_key, file_key="phonon_file"
        )

    @property
    def md_path(self) -> str | None:
        """MD prediction path, downloading when a URL is present."""
        return self._metric_pred_path("md")

    @property
    def diatomics_path(self) -> str | None:
        """Diatomics prediction path, downloading when a URL is present."""
        return self._metric_pred_path("diatomics")

    @property
    def is_active(self) -> bool:
        """Return whether the model participates in current benchmark views."""
        return self.metadata["lifecycle"] == "active"

    @classmethod
    def active(cls) -> tuple[Self, ...]:
        """Active models included by default in plots, metrics, and eval scripts."""
        return tuple(model for model in cls if model.is_active)

    @classmethod
    def from_ref(cls, ref: str | Self) -> Self:
        """Resolve an enum name, normalized value, or canonical model key."""
        if isinstance(ref, cls):
            return ref
        if member := cls.__members__.get(ref):
            return member
        for model in cls:
            if ref == model.key:
                return model
        if member := cls._missing_(ref):
            return member
        raise ValueError(f"{ref!r} not found in Model enum names or keys")

    @classmethod
    def _missing_(cls, value: object) -> Self | None:
        """Normalize casing and punctuation before matching enum values.

        Lets CLI arguments like --models mace-mp-0 resolve to mace_mp_0.
        """
        if isinstance(value, str):
            converted_value = re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")
            return cls._value2member_map_.get(converted_value)  # ty: ignore[invalid-return-type]

        return None


class DataFiles(Files):
    """Data-file registry whose paths and metadata come from ``data-files.yml``."""

    def __new__(cls, val: str) -> Self:
        """Create a registry member identified only by its stable YAML key."""
        obj = str.__new__(cls, val)
        obj._value_ = val
        return obj

    all_mp_tasks = auto()
    mp_computed_structure_entries = auto()
    mp_elemental_ref_entries = auto()
    mp_energies = auto()
    mp_patched_phase_diagram = auto()
    mp_trj_original = auto()
    mp_trj_extxyz = auto()
    wbm_computed_structure_entries = auto()
    wbm_relaxed_atoms = auto()
    wbm_initial_structures = auto()
    wbm_initial_atoms = auto()
    wbm_summary = auto()
    phonondb_pbe_103_structures = auto()
    phonondb_pbe_103_kappa_no_nac = auto()
    wbm_dft_geo_opt_symprec_1e_2 = auto()
    wbm_dft_geo_opt_symprec_1e_5 = auto()
    dynamat_v1_0_md_trajectories = auto()
    diatomics_dft_reference = auto()

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
    def rel_path(self) -> str:
        """Path from the canonical ``data-files.yml`` registry entry."""
        rel_path = self.yaml[self.name].get("path")
        if not isinstance(rel_path, str):
            raise TypeError(f"{self.name!r} does not have a string path")
        return rel_path

    @property
    def path(self) -> str:
        """Return the downloaded file only after validating its declared MD5."""
        abs_path = f"{type(self).base_dir}/{self.rel_path}"
        return resolve_verified_file(
            abs_path,
            url=self.url,
            label=self.name,
            md5=self.yaml[self.name].get("md5"),
        )
