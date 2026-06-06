// Shared ref-line idioms for the inline matterviz figures in site/src/routes.
import type { LegendConfig, RefLine } from 'matterviz/plot'

export const dashed = { dash: `4`, width: 1 }

// Long legends (many series) overlap the data when floated inside the plot, so render them
// as a full-width row of content-sized items below it. flex-wrap sizes each item to its
// label and wraps, avoiding the cross-column overlap a fixed-width grid causes. matterviz
// auto-places a wide legend in a reserved bottom margin.
export const wide_legend: LegendConfig = {
  draggable: false,
  filterable: false, // hide the search input; the legend is a static reference
  style: `width: 100%; display: flex; flex-wrap: wrap; justify-content: center; gap: 2px 14px;`,
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
