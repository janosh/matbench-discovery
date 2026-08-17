// Shared ref-line idioms for the inline matterviz figures in site/src/routes.
import { MODELS } from '$lib/models.svelte'
import type { LegendConfig, RefLine } from 'matterviz/plot'
import type { Attachment } from 'svelte/attachments'
import { SvelteSet } from 'svelte/reactivity'
import { dismiss_on_outside_press } from 'svelte-widgets/attachments'

// === multi-model payload styling ===
// .jsonl payloads hold position-independent data only (no color/order) so models merge
// cleanly. The json_payload plugin (vite.config.ts) runs each through attach_style on
// import, so pages get models pre-colored in discovery-F1-desc leaderboard order;
// order_models re-sorts the few figures wanting a different order (MAE, AUC, sigma).

// stable color + discovery metrics (unique-prototypes test set), keyed only by the
// immutable model_key. Display labels are presentation and never identify records.
const model_meta: Record<string, { color?: string; f1: number; mae: number }> = {}
for (const model of MODELS) {
  const metrics = model.metrics?.discovery?.unique_prototypes
  model_meta[model.model_key] = {
    color: model.color,
    f1: metrics?.F1 ?? -Infinity,
    mae: metrics?.MAE ?? Infinity,
  }
}
const meta = (model: { model_key: string }) => model_meta[model.model_key]
// discovery MAE key (ascending = best first) for the hull-dist box + rolling figures
export const model_mae = (model: { model_key: string }): number =>
  meta(model)?.mae ?? Infinity

// non-mutating re-sort by a numeric key (ascending = first), e.g. model_mae or
// (m) => -m.auc. color is attached by attach_style, so this only reorders.
export const order_models = <T>(models: T[], order: (model: T) => number): T[] =>
  models.toSorted((row_a, row_b) => order(row_a) - order(row_b))

// import-time pass (run by the json_payload plugin): attach each model's stable color
// and sort into discovery-F1-desc leaderboard order. Deriving it on import keeps the
// committed .jsonl position-independent and merge-friendly.
export const attach_style = <
  T extends { model_key: string; label: string },
  P extends { models: T[] },
>(
  payload: P,
) => ({
  ...payload,
  models: order_models(
    payload.models.map((model) => ({ ...model, color: meta(model)?.color })),
    (model) => -(meta(model)?.f1 ?? -Infinity),
  ),
})

export const dashed = { dash: `4`, width: 1 }

// Default series colors across data/model histograms
// (blue = primary/stable series, red = secondary/unstable series)
export const series_blue = `#636efa`
export const series_red = `#ef553b`

// Long legends (many series) overlap the data when floated inside the plot, so render them
// as a full-width row of content-sized items below it. flex-wrap sizes each item to its
// label and wraps, avoiding the cross-column overlap a fixed-width grid causes. matterviz
// auto-places a wide legend in a reserved bottom margin.
export const wide_legend: LegendConfig = {
  draggable: false,
  filterable: false, // hide the search input; the legend is a static reference
  style: `width: 100%; display: flex; flex-wrap: wrap; justify-content: center; gap: 2px 14px;`,
}

// Collapsible per-model legend: group headers collapse/expand while items still toggle
// individual model visibility. Attach `collapse_on_outside_click` to the plot wrapper.
export const make_models_legend = (legend_group = `Toggle Models`) => {
  const collapsed_groups = new SvelteSet([legend_group])
  const toggle_group = (group: string) => {
    if (!collapsed_groups.delete(group)) collapsed_groups.add(group)
  }
  const toggle = () => toggle_group(legend_group)
  const legend: LegendConfig = {
    ...wide_legend,
    collapsed_groups,
    // keep expanded legend readable over plot points
    style: `${wide_legend.style} --plot-legend-bg-color: light-dark(rgb(255, 255, 255), rgb(40, 40, 40))`,
    on_group_toggle: toggle_group,
  }
  // Only the legend counts as inside — a click on the plot or the controls above it
  // collapses too — so pass the surface as `inside` rather than attaching to `node`
  // (which click_outside would treat as inside). `scope` keeps a sibling figure's
  // legend from counting. `release` listens for click, not pointerdown, matching the
  // legend's own toggles, and the listener is capture-phase either way, so it still
  // sees clicks whose inner handlers stop propagation.
  const collapse_on_outside_click: Attachment = (node) =>
    dismiss_on_outside_press({
      inside: [`.legend`],
      scope: node,
      dismiss_on: `release`,
      callback: () => collapsed_groups.add(legend_group),
    })
  return { legend_group, legend, collapse_on_outside_click, toggle }
}

// full-span y=x parity diagonal; a diagonal ref line is clipped to the axis range
// (unlike a 2-point data series, it always reaches both plot corners)
export const parity_diagonal: RefLine = {
  type: `diagonal`,
  slope: 1,
  intercept: 0,
  label: `DFT = ML`,
  show_in_legend: true,
  style: { dash: `4 4`, color: `gray` },
}

// dashed vertical line whose label sits in the top margin, above the plot area
export const labeled_vline = (x: number, text: string): RefLine => ({
  type: `vertical`,
  x,
  style: dashed,
  annotation: { text, gap: 0, edge_padding: 0, offset: { y: -5 } },
})

// line-less label centered mid-plot (e.g. stable/unstable region labels)
export const floating_label = (x: number, text: string, color: string): RefLine => ({
  type: `vertical`,
  x,
  style: { width: 0 },
  annotation: { text, position: `center`, font_size: `15px`, color },
})
