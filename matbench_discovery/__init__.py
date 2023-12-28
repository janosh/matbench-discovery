"""Global variables used all across the matbench_discovery package."""

import json
import os
import warnings
from datetime import datetime

import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
from pymatviz.utils import styled_html_tag

ROOT = os.path.dirname(os.path.dirname(__file__))  # repo root directory
DATA_DIR = f"{ROOT}/data"  # directory to store raw data
SITE_FIGS = f"{ROOT}/site/src/figs"  # directory for interactive figures
SITE_MODELS = f"{ROOT}/site/src/routes/models"  # directory to write model analysis
FIGSHARE = f"{ROOT}/data/figshare"
SCRIPTS = f"{ROOT}/scripts"
PDF_FIGS = f"{ROOT}/paper/figs"  # directory for light-themed PDF figures

for directory in [SITE_FIGS, SITE_MODELS, FIGSHARE, PDF_FIGS]:
    os.makedirs(directory, exist_ok=True)

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

id_col = "material_id"
init_struct_col = "initial_structure"
struct_col = "structure"
e_form_col = "formation_energy_per_atom"
formula_col = "formula"
stress_col = "stress"
stress_trace_col = "stress_trace"
n_sites_col = "n_sites"
entry_col = "computed_structure_entry"

# load figshare 1.0.0
with open(f"{FIGSHARE}/1.0.0.json") as file:
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
