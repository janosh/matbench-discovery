<script lang="ts">
  import type { Snippet } from 'svelte'

  // snippets the page passes in, rendered below where each belongs
  let {
    task_navigation,
    min_relaxed_structures,
    geo_opt_metrics_table,
    model_comparison_scatter,
    diagnostic_model_picker,
    struct_rmsd_cdf_models,
    sym_ops_diff_bar,
    spg_sankeys,
    test_set,
  }: Record<string, Snippet | undefined> = $props()
</script>

# MLFF Geometry Optimization

This task measures how closely machine-learning force-field relaxations reproduce DFT-relaxed crystal structures and retain their symmetry.

Reference data: the [WBM test set](/benchmarks/geo-opt#test-set). The headline metric is RMSD AUC: the area under each model's cumulative RMSD distribution up to 0.05, with missing and unmatched structures scoring 0 (higher is better, 1 is perfect).

{@render task_navigation?.()}

## Leaderboard

Not all models relaxed every structure, but each reported model was evaluated on at least {@render min_relaxed_structures?.()} relaxations.

> **Symmetry-tolerance caveat:** RMSD is `symprec`-invariant, so RMSD AUC is reported once per model. The leaderboard shows symmetry metrics at both `symprec=1e-2` and `symprec=1e-5`, while Aggregate Diagnostics use `symprec=1e-5`. Compare symmetry values within a single view.

{@render geo_opt_metrics_table?.()}

## Model Comparison

Use the axis, color, and size controls to compare the models with geometry-optimization metrics directly from the current model YAML data.

{@render model_comparison_scatter?.()}

## Aggregate Diagnostics

These views use the per-structure `symprec=1e-5` analyses. The model picker applies to every view below.

{@render diagnostic_model_picker?.()}

### Cumulative Distribution of RMSD

{@render struct_rmsd_cdf_models?.()}

> Cumulative distribution of RMSD between ML- and DFT-relaxed structures, as a fraction of all WBM structures. A curve plateaus below 100% by the share of structures the model failed to relax or that can't be matched to DFT. The AUC in each legend is the model's RMSD AUC: the area under its curve up to 0.05, divided by 0.05.

### Difference in Number of Symmetry Operations vs DFT

{@render sym_ops_diff_bar?.()}

> Difference in number of symmetry operations of ML vs DFT-relaxed structures. Models are sorted by the standard deviation σ of ΔN<sub>sym ops</sub> = N<sub >sym ops,ML</sub> - N<sub>sym ops,DFT</sub>.

### Sankey Diagrams for ML vs DFT Spacegroups

The Sankey diagrams show corresponding spacegroups of DFT-relaxed and MLFF-relaxed structures at `symprec=1e-5`. For visual clarity, only the 10 most common pairs of (DFT, MLFF) spacegroups are shown.

{@render spg_sankeys?.()}

{@render test_set?.()}

## Methodology

**RMSD** is the normalized, unitless RMSD between ML- and DFT-relaxed ground-state structures returned by `pymatgen`'s `StructureMatcher` (`stol=1`, `ltol=0.2`, `angle_tol=5°`, `scale=False`). It measures atomic positions only: a strained cell within `ltol` and `angle_tol` matches with the same RMSD as an unstrained one, so volume and cell-shape errors below those tolerances go unpenalized. **RMSD AUC** averages max(1 - RMSD/0.05, 0) over all 256,963 WBM structures, which equals the area under the cumulative RMSD distribution up to 0.05, divided by 0.05. Missing predictions and structures whose cell shape or atomic arrangement differs from DFT's beyond the matcher's tolerances score 0, like an RMSD ≥ 0.05. Looser tolerances (`ltol=0.5`, `angle_tol=20°`) pair most of them, but at a median RMSD of 0.11–0.15: genuinely different structures, not near misses.

[RMSD AUC replaced mean RMSD in November 2026](https://github.com/janosh/matbench-discovery/pull/418). Mean RMSD set unmatched structures to 1.0, about 100× a typical matched RMSD, so failure penalties made up 71–81% of every model's value. It also ignored structures a model never predicted, and the old CDF plot dropped unmatched ones altogether. Metric, CDF plot and CPS now share one full-test-set definition, with RMSD AUC entering CPS directly instead of rescaled against a fixed 0.15 RMSD baseline.

The 0.05 cutoff sets what counts as a different geometry. `StructureMatcher` divides the RMS displacement by (V/n)<sup>1/3</sup>, which is 2.3–3.2 Å across WBM (5th–95th percentile, median 2.7 Å), so 0.05 is an RMS displacement of about 0.13 Å (0.11–0.16 Å), similar to room-temperature thermal vibrations in solids. The median matched RMSD is 0.0045 and only 7% of matched structures exceed 0.05. Rankings barely depend on the cutoff: Kendall τ between the ranking at 0.05 and those at 0.03 and 0.10 is ≥ 0.96. Some cutoff is unavoidable: at 1 (the matcher's `stol`), RMSD AUC reduces to 1 minus the old mean RMSD with missing structures also counted as 1.0, and cutoff-free alternatives such as the median RMSD ignore how often a model fails. About 3.7k structures can't be matched for any of the 49 active models because DFT relaxation changed them substantially (only 0.5% of their initial structures match their DFT-relaxed ones, against 90% across WBM). They lower every score without changing the ranking.

**Optimizer**, **Steps**, **f<sub>max</sub>**, and **Filter** show the ASE optimizer, maximum relaxation steps, force convergence criterion (eV/Å), and cell filter used during relaxation.

Symmetry detection uses [`moyopy`](https://github.com/spglib/moyo), a Rust successor to [`spglib`](https://spglib.readthedocs.io), with Hall numbers in `spglib`'s setting. Model analysis filenames record the `moyopy` version. Σ<sub>=</sub> denotes the fraction of structures with identical ML- and DFT-relaxed spacegroups. Σ<sub>↓</sub> / Σ<sub>↑</sub> denote the fractions whose numbers of symmetry operations decreased / increased during MLFF relaxation. Symmetry metrics depend on the `symprec` tolerance, so results are shown at multiple values. See the [`spglib` docs](https://spglib.readthedocs.io/en/latest/variable.html#symprec) and [paper](https://arxiv.org/html/1808.01590v2) for details.

The WBM DFT references were generated with `MPRelaxSet`: [`ISYM=2`](https://vasp.at/wiki/index.php/ISYM), [`ISIF=3`](https://vasp.at/wiki/index.php/ISIF), and `IBRION=2` (conjugate gradient). Most MLFF relaxations instead use `FIRE`. Different symmetry constraints and optimizers can reach different minima, so a lower symmetry-match rate can occasionally indicate a valid symmetry-broken structure rather than a model error. See [`MPRelaxSet.yaml`](https://github.com/materialsproject/pymatgen/blob/bf2cd24b647a33/src/pymatgen/io/vasp/MPRelaxSet.yaml#L10). Thanks to [Alex Ganose](https://scholar.google.co.uk/citations?user=nVJFXWwAAAAJ) for highlighting this distinction.

<style>
  h3 {
    text-align: center;
  }
</style>
