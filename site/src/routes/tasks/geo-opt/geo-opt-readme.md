# MLFF Geometry Optimization Analysis

> Disclaimer: There is a caveat to the structure similarity analysis below. The WBM test set was generated using the `MPRelaxSet` which uses [`ISYM=2`](https://vasp.at/wiki/index.php/ISYM) (VASP's default), [`ISIF=3`](https://vasp.at/wiki/index.php/ISIF) (see [`MPRelaxSet.yaml`](https://github.com/materialsproject/pymatgen/blob/bf2cd24b647a33/src/pymatgen/io/vasp/MPRelaxSet.yaml#L10)) and `IBRION=2` (i.e. conjugate gradient). `ISIF=3` means positions, cell shape and cell volume are all allowed to change. However, `ISYM=2` means the charge density is symmetrized according to the cell symmetry. Unsure how much this prevents the structure from breaking symmetry during relaxation. In any case, the optimizers are different. MLIP's mostly use `FIRE` which has a momentum term whereas VASP uses conjugate gradient to update ion positions. `FIRE`, on average, finds lower ground states than conjugate gradients which in some cases may require symmetry breaking to reach. This is not a mistake of the model and so higher Σ<sub>=</sub> (the percentage of structures with matching ML and DFT spacegroups) is not strictly speaking indicative of a better model (though for the purpose of this benchmark it probably is). Thanks to [Alex Ganose](https://scholar.google.co.uk/citations?user=nVJFXWwAAAAJ) for pointing this out! Undiscovered lower energy and lower symmetry structures in the well-explored chemical systems covered by WBM and MP are not expected to be a common occurrence. Hence we believe this analysis still provides some useful insight.

All plots/metrics below evaluate the quality of MLFF relaxations for the 257k crystal structures in the [WBM test set](https://nature.com/articles/s41524-020-00481-6). Not all models were able to relax all structures (user/cluster error may explain some failures) but every model was evaluated on at least <slot name="min_relaxed_structures" /> relaxations.

Symmetry detection was performed with the excellent Rust library [`moyopy`](https://github.com/spglib/moyo), a ~4x faster successor to [`spglib`](https://spglib.readthedocs.io).

<slot name="geo_opt_metrics_table" />

> **RMSD** measures the root-mean-square displacement between ML and DFT ground state structures. **Optimizer**, **Steps**, **f<sub>max</sub>**, and **Filter** show the ASE optimizer, maximum relaxation steps, force convergence criterion (eV/Å), and cell filter used during structure relaxation. Σ<sub>=</sub> / Σ<sub>↓</sub> / Σ<sub>↑</sub> denote the fraction of structures that retain, increase, or decrease the symmetry of the DFT-relaxed structure during MLFF relaxation. The match criterion is for the ML ground state to have identical spacegroup as DFT. For Σ<sub>↓</sub> / Σ<sub>↑</sub>, the number of symmetry operations for a structure increased / decreased during MLFF relaxation. Note that the symmetry metrics are sensitive to the `symprec` value passed to `spglib` so we show results for multiple values. See the [`spglib` docs](https://spglib.readthedocs.io/en/latest/variable.html#symprec) and [paper](https://arxiv.org/html/1808.01590v2) for details.

## Cumulative Distribution of RMSD

<slot name="struct_rmsd_cdf_models" />

> Cumulative distribution of RMSD between ML and DFT-relaxed structures.

## Difference in Number of Symmetry Operations vs DFT

<slot name="sym_ops_diff_bar" />

> Difference in number of symmetry operations of ML vs DFT-relaxed structures. Models are sorted by the standard deviation σ of ΔN<sub>sym ops</sub> = N<sub >sym ops,ML</sub> - N<sub>sym ops,DFT</sub>.

## Sankey Diagrams for ML vs DFT Spacegroups

The sankey diagrams show corresponding spacegroups of DFT-relaxed and MLFF-relaxed structures at `symprec=1e-5`. For visual clarity, only the 10 most common pairs of (DFT, MLFF) spacegroups are shown.
