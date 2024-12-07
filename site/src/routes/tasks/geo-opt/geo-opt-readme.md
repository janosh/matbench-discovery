# MLFF Geometry Optimization Analysis

> Disclaimer: There is a caveat to the structure similarity analysis  below. The WBM test set was generated using the `MPRelaxSet` which applies [`ISYM=2`](https://vasp.at/wiki/index.php/ISIF). This fixes the structure's symmetry. The MLFFs by contrast use the FIRE or LBFGS optimizers with no symmetry constraints. They may therefore in some cases relax to lower energy states with different symmetry. This is not a mistake of the model and so higher σ<sub>match</sub> (the percentage of structures with matching ML and DFT spacegroups) is not necessarily indicative of a better model. Thanks to Alex Ganose for pointing this out! Undiscovered lower energy structures in the relatively well-explored chemical systems covered by WBM and MP are not expected to be a common occurrence. Hence we believe this analysis still provides some signal so we left it on this secluded page.

All plots/metrics below evaluate the quality of MLFF relaxations for the 257k crystal structures in the [WBM test set](https://nature.com/articles/s41524-020-00481-6). Not all models were able to relax all structures (user/cluster error may explain some failures) but every model was evaluated on at least <slot name="min-relaxed-structures"/> relaxations.

Symmetry detection was performed with the excellent Rust library [`moyopy`](https://github.com/janosh/moyopy), a ~4x faster successor to the already outstanding [`spglib`](https://spglib.readthedocs.io).

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
