// Shared ref-line idioms for the inline matterviz figures in site/src/routes.
import type { LegendConfig, RefLine } from 'matterviz/plot'
import type { Attachment } from 'svelte/attachments'
import { SvelteSet } from 'svelte/reactivity'

export const dashed = { dash: `4`, width: 1 }

// plotly's first two qualitative colors, reused as the default series colors across the
// data/model histograms (blue = primary/stable series, red = secondary/unstable series)
export const plotly_blue = `#636efa`
export const plotly_red = `#ef553b`

// Long legends (many series) overlap the data when floated inside the plot, so render them
// as a full-width row of content-sized items below it. flex-wrap sizes each item to its
// label and wraps, avoiding the cross-column overlap a fixed-width grid causes. matterviz
// auto-places a wide legend in a reserved bottom margin.
export const wide_legend: LegendConfig = {
  draggable: false,
  filterable: false, // hide the search input; the legend is a static reference
  style: `width: 100%; display: flex; flex-wrap: wrap; justify-content: center; gap: 2px 14px;`,
}

type ModelsLegendConfig = {
  legend_group: string
  legend: LegendConfig
  collapse_on_outside_click: Attachment
}

// Collapsible per-model legend: group headers collapse/expand while items still toggle
// individual model visibility. Attach `collapse_on_outside_click` to the plot wrapper.
export const make_models_legend = (
  legend_group = `Toggle Models`,
): ModelsLegendConfig => {
  const collapsed_groups = new SvelteSet([legend_group])
  const legend: LegendConfig = {
    ...wide_legend,
    collapsed_groups,
    // keep expanded legend readable over plot points
    style: `${wide_legend.style} --plot-legend-bg-color: light-dark(rgb(255, 255, 255), rgb(40, 40, 40))`,
    on_group_toggle: (group) => {
      if (!collapsed_groups.delete(group)) collapsed_groups.add(group)
    },
  }
  // Capture-phase listener also sees clicks whose inner handlers stop propagation.
  const collapse_on_outside_click: Attachment = (node) => {
    const on_click = (event: MouseEvent) => {
      if (collapsed_groups.has(legend_group)) return
      const target = event.target
      if (!(target instanceof Node)) return
      const legend_el = node.querySelector(`.legend`)
      if (legend_el && !legend_el.contains(target)) {
        collapsed_groups.add(legend_group)
      }
    }
    document.addEventListener(`click`, on_click, true)
    return () => document.removeEventListener(`click`, on_click, true)
  }
  return { legend_group, legend, collapse_on_outside_click }
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
