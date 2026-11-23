# MLFF Geometry Optimization

This task measures how closely machine-learning force-field relaxations reproduce DFT-relaxed crystal structures and retain their symmetry.

Reference data: the [WBM test set](/benchmarks/geo-opt#test-set). The headline metric is RMSD AUC: the area under each model's cumulative RMSD distribution up to 0.05, with missing and unmatched structures scoring 0 (higher is better, 1 is perfect).

<slot name="task_navigation" />

## Leaderboard

Not all models relaxed every structure, but each reported model was evaluated on at least <slot name="min_relaxed_structures" /> relaxations.

> **Symmetry-tolerance caveat:** RMSD is `symprec`-invariant, so RMSD AUC is reported once per model. The leaderboard shows symmetry metrics at both `symprec=1e-2` and `symprec=1e-5`, while Aggregate Diagnostics use `symprec=1e-5`. Compare symmetry values within a single view.

<slot name="geo_opt_metrics_table" />

## Model Comparison

Use the axis, color, and size controls to compare the models with geometry-optimization metrics directly from the current model YAML data.

<slot name="model_comparison_scatter" />

## Aggregate Diagnostics

These views use the per-structure `symprec=1e-5` analyses. The model picker applies to every view below.

<slot name="diagnostic_model_picker" />

### Cumulative Distribution of RMSD

<slot name="struct_rmsd_cdf_models" />

> Cumulative distribution of RMSD between ML- and DFT-relaxed structures, as a fraction of all WBM structures. A curve plateaus below 100% by the share of structures the model failed to relax or that can't be matched to DFT. The AUC in each legend is the model's RMSD AUC: the area under its curve up to 0.05, divided by 0.05.

### Difference in Number of Symmetry Operations vs DFT

<slot name="sym_ops_diff_bar" />

> Difference in number of symmetry operations of ML vs DFT-relaxed structures. Models are sorted by the standard deviation σ of ΔN<sub>sym ops</sub> = N<sub >sym ops,ML</sub> - N<sub>sym ops,DFT</sub>.

### Sankey Diagrams for ML vs DFT Spacegroups

The Sankey diagrams show corresponding spacegroups of DFT-relaxed and MLFF-relaxed structures at `symprec=1e-5`. For visual clarity, only the 10 most common pairs of (DFT, MLFF) spacegroups are shown.

<slot name="spg_sankeys" />

<slot name="test_set" />

## Methodology

**RMSD** is the normalized, unitless structure-matching RMSD between ML- and DFT-relaxed ground-state structures, as returned by `pymatgen`'s `StructureMatcher` after matching (`stol=1`, `ltol=0.2`, `angle_tol=5°`, `scale=False`). It measures atomic positions only: a cell whose lattice differs from DFT's but stays within `ltol` and `angle_tol` matches with the same RMSD as an unstrained one, so volume and cell-shape errors below those tolerances are not penalized. **RMSD AUC** scores each structure as max(1 - RMSD/0.05, 0) and averages over all 256,963 WBM structures. That average equals the area under the cumulative RMSD distribution up to 0.05, divided by 0.05. A missing prediction, or a structure whose cell shape or atomic arrangement differs from DFT's beyond the matcher's tolerances, scores 0, the same as an RMSD of 0.05 or more. Loosening the tolerances (`ltol=0.5`, `angle_tol=20°`) pairs most such structures, but at a median RMSD of 0.11–0.15, so they are genuinely different structures rather than near misses.

RMSD AUC replaced mean RMSD in November 2026. The old mean assigned unmatched structures an RMSD of 1.0, about 100 times a typical matched RMSD. Those penalties made up 71–81% of every model's mean, so the metric mostly counted failures. It also ignored structures a model never predicted, and the old CDF plot dropped unmatched structures altogether. All three now use the same full-test-set definition, and RMSD AUC feeds CPS directly instead of being rescaled against a fixed 0.15 RMSD baseline. Model rankings barely depend on the 0.05 cutoff: Kendall τ between the rankings at cutoffs of 0.03 and 0.10 and the ranking at 0.05 is at least 0.96. About 3.7k structures can't be matched for any of the 49 active models: DFT relaxation changed them substantially (only 0.5% of their initial structures match their DFT-relaxed ones, against 90% across WBM). They lower every model's score but don't change the ranking. **Optimizer**, **Steps**, **f<sub>max</sub>**, and **Filter** show the ASE optimizer, maximum relaxation steps, force convergence criterion (eV/Å), and cell filter used during structure relaxation.

Symmetry detection uses [`moyopy`](https://github.com/spglib/moyo), a Rust successor to [`spglib`](https://spglib.readthedocs.io), with Hall numbers in `spglib`'s setting. Model analysis filenames record the `moyopy` version. Σ<sub>=</sub> denotes the fraction of structures with identical ML- and DFT-relaxed spacegroups. Σ<sub>↓</sub> / Σ<sub>↑</sub> denote the fractions whose numbers of symmetry operations decreased / increased during MLFF relaxation. Symmetry metrics depend on the `symprec` tolerance, so results are shown at multiple values. See the [`spglib` docs](https://spglib.readthedocs.io/en/latest/variable.html#symprec) and [paper](https://arxiv.org/html/1808.01590v2) for details.

The WBM DFT references were generated with `MPRelaxSet`: [`ISYM=2`](https://vasp.at/wiki/index.php/ISYM), [`ISIF=3`](https://vasp.at/wiki/index.php/ISIF), and `IBRION=2` (conjugate gradient). Most MLFF relaxations instead use `FIRE`. Different symmetry constraints and optimizers can reach different minima, so a lower symmetry-match rate can occasionally indicate a valid symmetry-broken structure rather than a model error. See [`MPRelaxSet.yaml`](https://github.com/materialsproject/pymatgen/blob/bf2cd24b647a33/src/pymatgen/io/vasp/MPRelaxSet.yaml#L10). Thanks to [Alex Ganose](https://scholar.google.co.uk/citations?user=nVJFXWwAAAAJ) for highlighting this distinction.

<style>
  h3 {
    text-align: center;
  }
</style>
