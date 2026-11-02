# GXTNet

GXTNet is a GPS-style graph transformer interatomic potential developed at LG AI Research: per layer a
GINE message-passing branch with edge updates and a sparse, distance-biased attention branch over the
same 6 A radius graph plus one virtual node per structure. Energies are conservative (forces and stress by
autograd), with a ZBL short-range repulsion.

| model | training set | checkpoint | notes |
| --- | --- | --- | --- |
| `gxtnet-s-mptrj` | MPtrj | 7,414,273 parameters | virtual node without centroid distance, attention cutoff envelope, dropout-free |

Model code and weights are not released; this is a predictions-only submission. The artifacts were produced
by the authors with the shared Matbench Discovery runners (float64 eager, TF32 off; diatomics float32) and an
identical-protocol MD driver. Metrics and cost provenance are in `gxtnet-s-mptrj.yml`. All prediction artifacts
are archived on Figshare under CC BY 4.0: https://doi.org/10.6084/m9.figshare.34885332.
