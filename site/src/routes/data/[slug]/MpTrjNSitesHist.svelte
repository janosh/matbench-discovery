<script lang="ts">
  // MPtrj sites-per-structure histogram with cumulative-share line on a secondary
  // axis; shared by the /data overview page and the MPtrj dataset page
  import mp_trj_hists from '$figs/mp-trj-hists.json.gz'
  import { dashed, series_blue } from '$lib/fig-helpers'
  import { BarPlot } from 'matterviz/plot'

  let { style = `height: 320px` }: { style?: string } = $props()

  const n_sites = mp_trj_hists[`n-sites`]
  // index where the cumulative share of structures crosses 90%
  const n_sites_90th = n_sites.x[n_sites.cumulative.findIndex((share) => share >= 0.9)]
</script>

<BarPlot
  series={[
    { ...n_sites, label: `Number of Structures`, color: series_blue },
    {
      x: n_sites.x,
      y: n_sites.cumulative,
      label: `Cumulative (%)`,
      render_mode: `line`,
      y_axis: `y2`,
    },
  ]}
  x_axis={{ range: [0, null] }}
  y_axis={{ label: `Number of Structures`, format: `~s` }}
  y2_axis={{ format: `.0%` }}
  ref_lines={[
    { type: `vertical`, x: n_sites_90th, style: dashed },
    { type: `horizontal`, y: 0.9, y_axis: `y2`, style: dashed },
  ]}
  show_legend
  legend={{ style: `left: auto; top: auto; right: 58px; bottom: 92px` }}
  show_controls={false}
  {style}
/>
