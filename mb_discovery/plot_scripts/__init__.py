import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio


px.defaults.labels = {
    "n_atoms": "Atom Count",
    "n_elems": "Element Count",
    "crystal_sys": "Crystal system",
    "spg_num": "Space group",
    "n_wyckoff": "Number of Wyckoff positions",
    "n_sites": "Lattice site count",
    "energy_per_atom": "Energy (eV/atom)",
    "e_form": "Formation energy (eV/atom)",
    "e_above_hull": "Energy above convex hull (eV/atom)",
    "e_above_hull_pred": "Predicted energy above convex hull (eV/atom)",
    "e_above_mp_hull": "Energy above MP convex hull (eV/atom)",
    "e_above_hull_error": "Error in energy above convex hull (eV/atom)",
}

pio.templates.default = "plotly_white"

# https://github.com/plotly/Kaleido/issues/122#issuecomment-994906924
# when seeing MathJax "loading" message in exported PDFs, try:
# pio.kaleido.scope.mathjax = None


plt.rc("font", size=16)
plt.rc("savefig", bbox="tight", dpi=200)
plt.rc("figure", dpi=200, titlesize=18)
plt.rcParams["figure.constrained_layout.use"] = True
