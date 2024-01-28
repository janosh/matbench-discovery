"""Global variables used all across the matbench_discovery package."""

import json
import os
import warnings
from datetime import datetime
from enum import StrEnum, unique
from importlib.metadata import Distribution

import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
from pymatviz.utils import styled_html_tag

pkg_name = "matbench-discovery"
direct_url = Distribution.from_name(pkg_name).read_text("direct_url.json") or "{}"
pkg_is_editable = json.loads(direct_url).get("dir_info", {}).get("editable", False)

PKG_DIR = os.path.dirname(__file__)
# repo root directory if editable install, else the pkg directory
ROOT = os.path.dirname(PKG_DIR) if pkg_is_editable else PKG_DIR
DATA_DIR = f"{ROOT}/data"  # directory to store raw data
SITE_FIGS = f"{ROOT}/site/src/figs"  # directory for interactive figures
# directory to write model analysis for website
SITE_MODELS = f"{ROOT}/site/src/routes/models"
SCRIPTS = f"{ROOT}/scripts"  # model and date analysis scripts
PDF_FIGS = f"{ROOT}/paper/figs"  # directory for light-themed PDF figures
FIGSHARE_DIR = f"{PKG_DIR}/figshare"

for directory in (SITE_FIGS, SITE_MODELS, PDF_FIGS):
    os.makedirs(directory, exist_ok=True)

os.makedirs(MP_DIR := f"{DATA_DIR}/mp", exist_ok=True)
os.makedirs(WBM_DIR := f"{DATA_DIR}/wbm", exist_ok=True)
# JSON files with URLS to data files on figshare

# directory to store model checkpoints downloaded from wandb cloud storage
CHECKPOINT_DIR = f"{ROOT}/wandb/checkpoints"
# wandb <entity>/<project name> to record new runs to
WANDB_PATH = "janosh/matbench-discovery"

# threshold on hull distance for a material to be considered stable
STABILITY_THRESHOLD = 0

timestamp = f"{datetime.now():%Y-%m-%d@%H-%M-%S}"
today = timestamp.split("@")[0]

# filter pymatgen warnings that spam the logs when e.g. applying corrections to
# ComputedStructureEntries or using PatchedPhaseDiagram to get e_above_hull
# warnings are:
# > No electronegativity for Ne. Setting to NaN. This has no physical meaning
# and MaterialsProject2020Compatibility to get formation energies
# > Failed to guess oxidation states for Entry
warnings.filterwarnings(action="ignore", category=UserWarning, module="pymatgen")


@unique
class Key(StrEnum):
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
class Task(StrEnum):
    """Thermodynamic stability prediction task types."""

    IS2RE = "IS2RE"  # initial structure to relaxed energy
    RS2RE = "RS2RE"  # relaxed structure to relaxed energy
    S2EFSM = "S2EFSM"  # structure to energy force stress magmom
    S2EFS = "S2EFS"  # structure to energy force stress
    S2RE = "S2RE"  # structure to relaxed energy (for models that learned a discrete
    # version of the PES like CGCNN+P)
    RP2RE = "RP2RE"  # relaxed prototype to relaxed energy
    IP2RE = "IP2RE"  # initial prototype to relaxed energy
    IS2E = "IS2E"  # initial structure to energy
    IS2RE_SR = "IS2RE-SR"  # initial structure to relaxed energy after ML relaxation


@unique
class Targets(StrEnum):
    """Thermodynamic stability prediction task types."""

    E = "E"
    EFS = "EFS"
    EFSM = "EFSM"


@unique
class ModelType(StrEnum):
    """Model types."""

    GNN = "GNN"
    UIP = "UIP-GNN"
    BO_GNN = "BO-GNN"
    Fingerprint = "Fingerprint"
    Transformer = "Transformer"


with open(f"{FIGSHARE_DIR}/1.0.0.json") as file:
    FIGSHARE_URLS = json.load(file)

# --- start global plot settings
ev_per_atom = styled_html_tag(
    "(eV/atom)", tag="span", style="font-size: 0.8em; font-weight: lighter;"
)
quantity_labels = dict(
    n_atoms="Atom Count",
    n_elems="Element Count",
    crystal_sys="Crystal system",
    spg_num="Space group",
    n_wyckoff="Number of Wyckoff positions",
    n_sites="Number of atoms",
    energy_per_atom=f"Energy {ev_per_atom}",
    e_form=f"DFT E<sub>form</sub> {ev_per_atom}",
    e_above_hull=f"E<sub>hull dist</sub> {ev_per_atom}",
    e_above_hull_mp2020_corrected_ppd_mp=f"DFT E<sub>hull dist</sub> {ev_per_atom}",
    e_above_hull_pred=f"Predicted E<sub>hull dist</sub> {ev_per_atom}",
    e_above_hull_mp=f"E<sub>above MP hull</sub> {ev_per_atom}",
    e_above_hull_error=f"Error in E<sub>hull dist</sub> {ev_per_atom}",
    vol_diff="Volume difference (A^3)",
    e_form_per_atom_mp2020_corrected=f"DFT E<sub>form</sub> {ev_per_atom}",
    e_form_per_atom_pred=f"Predicted E<sub>form</sub> {ev_per_atom}",
    material_id="Material ID",
    band_gap="Band gap (eV)",
    formula="Formula",
    stress="σ (eV/Å³)",  # noqa: RUF001
    stress_trace="1/3 Tr(σ) (eV/Å³)",  # noqa: RUF001
)
model_labels = dict(
    alignn="ALIGNN",
    alignn_ff="ALIGNN FF",
    alignn_pretrained="ALIGNN Pretrained",
    bowsr_megnet="BOWSR",
    chgnet="CHGNet",
    chgnet_megnet="CHGNet→MEGNet",
    cgcnn_p="CGCNN+P",
    cgcnn="CGCNN",
    m3gnet_megnet="M3GNet→MEGNet",
    m3gnet="M3GNet",
    m3gnet_direct="M3GNet DIRECT",
    m3gnet_ms="M3GNet MS",
    mace="MACE",
    megnet="MEGNet",
    megnet_rs2re="MEGNet RS2RE",
    voronoi_rf="Voronoi RF",
    wrenformer="Wrenformer",
    pfp="PFP",
    dft="DFT",
    wbm="WBM",
)
px.defaults.labels = quantity_labels | model_labels


global_layout = dict(
    paper_bgcolor="rgba(0,0,0,0)",
    font_size=13,
    # increase legend marker size and make background transparent
    legend=dict(itemsizing="constant", bgcolor="rgba(0, 0, 0, 0)"),
)
pio.templates["global"] = dict(layout=global_layout)
pio.templates.default = "plotly_dark+global"
px.defaults.template = "plotly_dark+global"

# https://github.com/plotly/Kaleido/issues/122#issuecomment-994906924
# when seeing MathJax "loading" message in exported PDFs,
# use pio.kaleido.scope.mathjax = None


plt.rc("font", size=14)
plt.rc("legend", fontsize=16, title_fontsize=16)
plt.rc("axes", titlesize=16, labelsize=16)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rc("figure", dpi=200, titlesize=16)
plt.rcParams["figure.constrained_layout.use"] = True
# --- end global plot settings
