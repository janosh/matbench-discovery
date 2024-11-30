# MLFF Geometry Optimization Analysis

All plots/metrics below evaluate the quality of MLFF relaxations for the 257k crystal structures in the [WBM test set](https://nature.com/articles/s41524-020-00481-6). Not all models were able to relax all structures (user/cluster error may explain some failures) but every model was evaluated on at least <slot name="min-relaxed-structures"/> relaxations.

<slot name="geo-opt-metrics-table"/>

> σ<sub>match</sub> / σ<sub>dec</sub> / σ<sub>inc</sub> denote the fraction of structures that retain, increase, or decrease the symmetry of the DFT-relaxed structure during MLFF relaxation. The match criterion is for the ML ground state to have identical spacegroup as DFT. For σ<sub>dec</sub> / σ<sub>inc</sub>, ML relaxation increased / decreased the set of symmetry operations on a structure. Note that the symmetry metrics are sensitive to the `symprec` value passed to `spglib` so we show results for multiple values. See the [`spglib` docs](https://spglib.readthedocs.io/en/latest/variable.html#symprec) and [paper](https://arxiv.org/html/1808.01590v2) for details.

<hr />

<slot name="struct-rmsd-cdf-models"/>

> Cumulative distribution of RMSD between ML and DFT-relaxed structures.

<hr />

<slot name="sym-ops-diff-bar"/>

> Difference in number of symmetry operations of ML vs DFT-relaxed structures. Models are sorted by the standard deviation σ of ΔN<sub>sym ops</sub> = N<sub >sym ops,ML</sub> - N<sub>sym ops,DFT</sub>.

---

## Sankey Diagrams for ML vs DFT Spacegroups

The sankey diagrams show corresponding spacegroups of DFT-relaxed and MLFF-relaxed structures. For visual clarity, only the 10 most common pairs of (DFT, MLFF) spacegroups are shown.
