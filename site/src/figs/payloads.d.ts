// Types for the data-only figure payloads in this directory. One exact ambient
// declaration per payload — these take precedence over the generic '*.json.gz' / '*.jsonl'
// fallbacks in app.d.ts, so each `import data from '$figs/<name>...'` is fully typed.
//
// Static payloads are committed as gzipped `<name>.json.gz`; multi-model payloads as
// line-delimited `<name>.jsonl` (one model per line, reassembled into the aggregate
// shape) so concurrent model submissions git-merge cleanly. Both are loaded at build by
// the json_payload plugin (vite.config.ts) and import as a parsed default export. The
// plugin pipes each .jsonl payload through attach_style, so models arrive in leaderboard
// (discovery-F1-desc) order with the `color` it injects, typed on KeyedModel and on the
// label-only payloads that read it.
//
// This file must stay import/export-free at the top level so the helper interfaces below
// are global and the module declarations stay ambient.

interface XY<TX = number> {
  x: TX[]
  y: number[]
}

// a per-model series: x/y arrays plus the model's display label
interface LabeledXY<TX = number> extends XY<TX> {
  label: string
}

// pre-binned histogram: x = bin centers, y = counts
interface HistBins extends XY {
  bar_width: number
}

// per-model entry in the discovery/diagnostic payloads: `key` (= MODELS model_key, used
// e.g. for the compliance join in discovery-metric-figs.md), a display `label`, and the
// stable `color` attach_style injects on import (undefined if the model isn't in MODELS;
// presentation, not committed data). Label-only payloads below inline color where read.
interface KeyedModel {
  key: string
  label: string
  color: string | undefined
}

// === tasks/discovery/tmi metrics ===
declare module '$figs/box-hull-dist-errors.jsonl' {
  const data: {
    // quantiles = [q05, q25, median, q75, q95] of each model's hull distance error
    models: (KeyedModel & { quantiles: number[] })[]
  }
  export default data
}

declare module '$figs/cumulative-precision-recall.jsonl' {
  const data: {
    n_stable: number // number of stable materials in the WBM test set
    models: (KeyedModel & {
      x: number[] // number of materials validated, ranked most to least stable
      precision: number[]
      recall: number[]
      // each model's end of ranking: [n predicted stable, precision, recall] there
      end: [number, number, number]
    })[]
  }
  export default data
}

declare module '$figs/roc-models.jsonl' {
  const data: {
    models: (KeyedModel & {
      auc: number
      fpr: number[]
      tpr: number[]
    })[]
  }
  export default data
}

declare module '$figs/rolling-mae-vs-hull-dist.jsonl' {
  const data: {
    x: number[] // shared E above hull values (eV/atom)
    models: (KeyedModel & {
      y: number[]
    })[]
    density: XY // rolling count of test-set structures per hull-dist bin (on y2)
  }
  export default data
}

declare module '$figs/hist-clf-pred-hull-dist.jsonl' {
  const data: {
    bin_centers: number[] // shared hull-dist bins (eV/atom)
    // per-model stability-classification counts per bin
    models: (KeyedModel & { f1: number } & Record<`tp` | `fn` | `fp` | `tn`, number[]>)[]
  }
  export default data
}

// === tasks/discovery/tmi extras ===
declare module '$figs/element-prevalence-vs-error.jsonl' {
  const data: {
    elements: string[] // element symbols, same order as occurrences
    occurrences: (number | null)[] // MP training-set occurrence count per element
    // mean error per element (color read by the per-element scatter)
    models: { label: string; color: string | undefined; y: (number | null)[] }[]
  }
  export default data
}

declare module '$figs/scatter-largest-fp-diff-each-error.jsonl' {
  const data: {
    fp_diff: number[] // shared |SSFP_initial - SSFP_final| values
    models: {
      label: string
      color: string | undefined
      mae: number
      y: (number | null)[]
    }[]
  }
  export default data
}

declare module '$figs/scatter-largest-each-errors-fp-diff.jsonl' {
  const data: { models: (LabeledXY & { mae: number })[] }
  export default data
}

declare module '$figs/hist-largest-each-errors-fp-diff.jsonl' {
  const data: {
    // fingerprint-diff histograms for each model's 100 worst (err_max) and best
    // (err_min) hull-dist predictions
    models: { label: string; err_min: HistBins; err_max: HistBins }[]
  }
  export default data
}

// === data pages ===
declare module '$figs/hist-wbm-e-form-per-atom.json.gz' {
  const data: HistBins
  export default data
}

declare module '$figs/hist-wbm-hull-dist.json.gz' {
  const data: {
    bar_width: number
    stable: XY
    unstable: XY
    mean: number // mean WBM hull distance (eV/atom)
    std: number
  }
  export default data
}

declare module '$figs/spacegroup-sunbursts.json.gz' {
  // flat plotly sunburst arrays; matterviz sunburst_from_labels_parents nests them
  interface SunburstArrays {
    labels: string[]
    parents: string[]
    values: number[]
    ids?: string[] // sunburst_data only emits ids when the plotly trace has them
  }
  const data: { mp: SunburstArrays; wbm: SunburstArrays }
  export default data
}

declare module '$figs/mp-vs-mp-trj-vs-wbm-arity-hist.json.gz' {
  // fraction of structures per number of elements in formula, by dataset
  const data: { datasets: (LabeledXY & { color: string })[] }
  export default data
}

declare module '$figs/mp-trj-hists.json.gz' {
  const data: {
    'e-form': HistBins
    forces: HistBins
    stresses: HistBins
    magmoms: HistBins
    'n-sites': HistBins & { cumulative: number[] } // cumulative share of structures
  }
  export default data
}

declare module '$figs/mp-elemental-ref-energies.json.gz' {
  // x = atomic number, y = lowest energy of any unary structure for that element
  const data: XY
  export default data
}

declare module '$figs/element-counts-mp-vs-wbm.json.gz' {
  // x = element symbols sorted by count, one series per dataset (WBM, MP)
  const data: {
    raw: (LabeledXY<string> & { color: string })[]
    normalized: (LabeledXY<string> & { color: string })[]
  }
  export default data
}

// === phonons ===
declare module '$figs/kappa-103-analysis.jsonl' {
  // per-material kappa-103 diagnostics vs the phononDB-PBE reference. All
  // per-material arrays are aligned to material_ids; null = material missing from
  // the model's predictions (or the value couldn't be computed)
  const data: {
    material_ids: string[]
    formulas: string[]
    spg_nums: number[] // spacegroup numbers (client derives crystal system)
    kappa_dft: (number | null)[] // DFT kappa_L at 300K (W/mK)
    models: (KeyedModel & {
      kappa_ml: (number | null)[] // ML kappa_L at 300K (W/mK)
      srme: (number | null)[] // per-material SRME in [0, 2]
      imag_modes: (boolean | null)[] // imaginary phonon modes after relaxation
      broken_sym: (boolean | null)[] // symmetry broken during relaxation
      max_steps: (boolean | null)[] // relaxation hit max steps (non-converged)
      freq_w1: (number | null)[] // phonon spectrum Wasserstein-1 dist vs DFT (THz)
      freq_w1_mean: number | null
      // ML vs DFT phonon frequency quantile-quantile parity pairs (THz), 17
      // quantile levels per material concatenated over all comparable materials
      freq_pairs: { dft: number[]; ml: number[] }
    })[]
  }
  export default data
}

// === geo-opt ===
declare module '$figs/struct-rmsd-cdf.jsonl' {
  const data: { models: (LabeledXY & { auc: number })[] }
  export default data
}

declare module '$figs/sym-ops-diff-bar.jsonl' {
  // histogram of symmetry-operation count changes during relaxation (symprec=1e-5)
  const data: { models: (LabeledXY & { sigma: number })[] }
  export default data
}

declare module '$figs/spg-sankeys.jsonl' {
  // DFT vs model spacegroup flows (symprec=1e-5); key matches MODELS model_key. flat
  // arrays; matterviz sankey_from_links(source, target, value, labels) builds the graph
  const data: {
    models: (KeyedModel & {
      labels: string[]
      source: number[]
      target: number[]
      value: number[]
    })[]
  }
  export default data
}
