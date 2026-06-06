// Re-export barrel for the data-only figure payloads in this directory, so pages can
// `import { roc_models } from '$figs'`. Payload shapes are declared per file in
// payloads.d.ts. These must stay PURE re-exports: a module-level binding (e.g. a cast)
// would tie all payloads into one chunk, whereas binding-level re-exports let the
// bundler put each payload only into the route chunks that actually use it.
export { default as box_hull_dist_errors } from '$figs/box-hull-dist-errors.json.gz'
export { default as cumulative_precision_recall } from '$figs/cumulative-precision-recall.json.gz'
export { default as element_counts_mp_vs_wbm } from '$figs/element-counts-mp-vs-wbm.json.gz'
export { default as element_prevalence_vs_error } from '$figs/element-prevalence-vs-error.json.gz'
export { default as hist_clf_pred_hull_dist } from '$figs/hist-clf-pred-hull-dist.json.gz'
export { default as hist_largest_each_errors_fp_diff } from '$figs/hist-largest-each-errors-fp-diff.json.gz'
export { default as hist_wbm_e_form_per_atom } from '$figs/hist-wbm-e-form-per-atom.json.gz'
export { default as hist_wbm_hull_dist } from '$figs/hist-wbm-hull-dist.json.gz'
export { default as mp_elemental_ref_energies } from '$figs/mp-elemental-ref-energies.json.gz'
export { default as mp_trj_hists } from '$figs/mp-trj-hists.json.gz'
export { default as arity_hist } from '$figs/mp-vs-mp-trj-vs-wbm-arity-hist.json.gz'
export { default as roc_models } from '$figs/roc-models.json.gz'
export { default as rolling_mae_vs_hull_dist } from '$figs/rolling-mae-vs-hull-dist.json.gz'
export { default as scatter_largest_each_errors_fp_diff } from '$figs/scatter-largest-each-errors-fp-diff.json.gz'
export { default as scatter_largest_fp_diff_each_error } from '$figs/scatter-largest-fp-diff-each-error.json.gz'
export { default as spacegroup_sunbursts } from '$figs/spacegroup-sunbursts.json.gz'
export { default as spg_sankeys } from '$figs/spg-sankeys.json.gz'
export { default as struct_rmsd_cdf } from '$figs/struct-rmsd-cdf.json.gz'
export { default as sym_ops_diff_bar } from '$figs/sym-ops-diff-bar.json.gz'
