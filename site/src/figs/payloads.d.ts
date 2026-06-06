// Types for the data-only figure payloads in this directory (written by
// matbench_discovery analysis scripts as gzipped JSON, decompressed at build time by
// the json_gz plugin in vite.config.ts). One exact ambient declaration per payload —
// these take precedence over the generic '*.json.gz' fallback in app.d.ts, so each
// `import data from '$figs/<name>.json.gz'` is fully typed with zero runtime
// indirection. This file must stay import/export-free at the top level so the
// helper interfaces below are global and the module declarations stay ambient.

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

// === models/tmi discovery metrics ===
declare module '$figs/box-hull-dist-errors.json.gz' {
  const data: {
    // quantiles = [q05, q25, median, q75, q95] of each model's hull distance error
    models: { label: string; color: string; quantiles: number[] }[]
  }
  export default data
}

declare module '$figs/cumulative-precision-recall.json.gz' {
  const data: {
    n_stable: number // number of stable materials in the WBM test set
    models: {
      label: string
      color: string
      x: number[] // number of materials validated, ranked most to least stable
      precision: number[]
      recall: number[]
      // each model's end of ranking: [n predicted stable, precision, recall] there
      end: [number, number, number]
    }[]
  }
  export default data
}

declare module '$figs/roc-models.json.gz' {
  const data: { models: { label: string; auc: number; fpr: number[]; tpr: number[] }[] }
  export default data
}

declare module '$figs/rolling-mae-vs-hull-dist.json.gz' {
  const data: {
    x: number[] // shared E above hull values (eV/atom)
    models: { label: string; color: string; y: number[]; visible?: boolean }[]
    density: XY // rolling count of test-set structures per hull-dist bin (on y2)
  }
  export default data
}

declare module '$figs/hist-clf-pred-hull-dist.json.gz' {
  const data: {
    bin_centers: number[] // shared hull-dist bins (eV/atom)
    // per-model stability-classification counts per bin
    models: ({ label: string; f1: number } & Record<
      `tp` | `fn` | `fp` | `tn`,
      number[]
    >)[]
  }
  export default data
}

// === models/tmi extras ===
declare module '$figs/element-prevalence-vs-error.json.gz' {
  const data: {
    elements: string[] // element symbols, same order as occurrences
    occurrences: number[] // MP training-set occurrence count per element
    models: { label: string; color: string; y: number[] }[] // mean error per element
  }
  export default data
}

declare module '$figs/scatter-largest-fp-diff-each-error.json.gz' {
  const data: {
    fp_diff: number[] // shared |SSFP_initial - SSFP_final| values
    models: { label: string; mae: number; color: string; y: number[] }[]
  }
  export default data
}

declare module '$figs/scatter-largest-each-errors-fp-diff.json.gz' {
  const data: { models: (LabeledXY & { mae: number })[] }
  export default data
}

declare module '$figs/hist-largest-each-errors-fp-diff.json.gz' {
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

// === geo-opt ===
declare module '$figs/struct-rmsd-cdf.json.gz' {
  const data: { models: (LabeledXY & { auc: number })[] }
  export default data
}

declare module '$figs/sym-ops-diff-bar.json.gz' {
  // histogram of symmetry-operation count changes during relaxation (symprec=1e-5)
  const data: { models: (LabeledXY & { sigma: number })[] }
  export default data
}

declare module '$figs/spg-sankeys.json.gz' {
  // DFT vs model spacegroup flows (symprec=1e-5); key matches MODELS model_key. flat
  // arrays; matterviz sankey_from_links(source, target, value, labels) builds the graph
  const data: {
    models: {
      key: string
      label: string
      labels: string[]
      source: number[]
      target: number[]
      value: number[]
    }[]
  }
  export default data
}
