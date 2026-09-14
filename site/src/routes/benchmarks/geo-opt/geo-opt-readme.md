# MLFF Geometry Optimization

This task measures how closely machine-learning force-field relaxations reproduce DFT-relaxed crystal structures and retain their symmetry.

Reference data: the [WBM test set](/benchmarks/geo-opt#test-set). The headline metric is normalized structure-matching RMSD (unitless, lower is better).

<slot name="task_navigation" />

## Leaderboard

Not all models relaxed every structure, but each reported model was evaluated on at least <slot name="min_relaxed_structures" /> relaxations.

> **Symmetry-tolerance caveat:** RMSD is `symprec`-invariant. The leaderboard shows symmetry metrics at both `symprec=1e-2` and `symprec=1e-5`, while Aggregate Diagnostics use `symprec=1e-5`. Compare symmetry values within a single view.

<slot name="geo_opt_metrics_table" />

## Model Comparison

Use the axis, color, and size controls to compare the models with geometry-optimization metrics directly from the current model YAML data.

<slot name="model_comparison_scatter" />

## Aggregate Diagnostics

These views use the per-structure `symprec=1e-5` analyses. The model picker applies to every view below.

<slot name="diagnostic_model_picker" />

### Cumulative Distribution of RMSD

<slot name="struct_rmsd_cdf_models" />

> Cumulative distribution of RMSD between ML and DFT-relaxed structures.

### Difference in Number of Symmetry Operations vs DFT

<slot name="sym_ops_diff_bar" />

> Difference in number of symmetry operations of ML vs DFT-relaxed structures. Models are sorted by the standard deviation σ of ΔN<sub>sym ops</sub> = N<sub >sym ops,ML</sub> - N<sub>sym ops,DFT</sub>.

### Sankey Diagrams for ML vs DFT Spacegroups

The Sankey diagrams show corresponding spacegroups of DFT-relaxed and MLFF-relaxed structures at `symprec=1e-5`. For visual clarity, only the 10 most common pairs of (DFT, MLFF) spacegroups are shown.

<slot name="spg_sankeys" />

<slot name="test_set" />

## Methodology

<details>
<summary>Structure matching, symmetry, and relaxation protocol</summary>

**RMSD** measures the normalized, unitless structure-matching RMSD between ML- and DFT-relaxed ground state structures, as returned by `pymatgen`'s `StructureMatcher` after matching. **Optimizer**, **Steps**, **f<sub>max</sub>**, and **Filter** show the ASE optimizer, maximum relaxation steps, force convergence criterion (eV/Å), and cell filter used during structure relaxation.

Symmetry detection uses [`moyopy`](https://github.com/spglib/moyo), a Rust successor to [`spglib`](https://spglib.readthedocs.io). Σ<sub>=</sub> denotes the fraction of structures with identical ML- and DFT-relaxed spacegroups. Σ<sub>↓</sub> / Σ<sub>↑</sub> denote the fractions whose numbers of symmetry operations decreased / increased during MLFF relaxation. Symmetry metrics depend on the `symprec` tolerance, so results are shown at multiple values. See the [`spglib` docs](https://spglib.readthedocs.io/en/latest/variable.html#symprec) and [paper](https://arxiv.org/html/1808.01590v2) for details.

The WBM DFT references were generated with `MPRelaxSet`: [`ISYM=2`](https://vasp.at/wiki/index.php/ISYM), [`ISIF=3`](https://vasp.at/wiki/index.php/ISIF), and `IBRION=2` (conjugate gradient). Most MLFF relaxations instead use `FIRE`. Different symmetry constraints and optimizers can reach different minima, so a lower symmetry-match rate can occasionally indicate a valid symmetry-broken structure rather than a model error. See [`MPRelaxSet.yaml`](https://github.com/materialsproject/pymatgen/blob/bf2cd24b647a33/src/pymatgen/io/vasp/MPRelaxSet.yaml#L10). Thanks to [Alex Ganose](https://scholar.google.co.uk/citations?user=nVJFXWwAAAAJ) for highlighting this distinction.

</details>

<style>
  h3 {
    text-align: center;
  }
</style>
