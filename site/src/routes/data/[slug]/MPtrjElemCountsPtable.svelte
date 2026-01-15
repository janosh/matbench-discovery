<script lang="ts">
  import { data_files, PtableHeatmap } from '$lib'
  import type { ElementSymbol } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'

  let {
    count_mode = `occurrence`,
    log = false,
    color_scale = `interpolateViridis`,
  }: {
    count_mode?: string
    log?: boolean
    color_scale?: D3InterpolateName
  } = $props()

  const elem_counts = import.meta.glob(`../mp-trj-element-counts-by-*.json`, {
    eager: true,
    import: `default`,
  }) as Record<string, Record<ElementSymbol, number>>

  let mp_trj_elem_counts = $derived(
    elem_counts[`../mp-trj-element-counts-by-${count_mode}.json`],
  )

  type DataFile = {
    figshare?: string
    url: string
    path: string
    description: string
  }
  const mp_trj_data = data_files[`mp_trj_json_gz`] as DataFile | undefined
  if (!mp_trj_data) {
    throw `mp_trj_json_gz not found in data-files.yml`
  }
</script>

<p>
  Element counts for
  <a href={mp_trj_data.figshare}>MPtrj training set</a> consisting of 1,580,395 structures
  which are frames of the DFT relaxations performed on all 154,719 MP materials.
</p>

<PtableHeatmap
  heatmap_values={mp_trj_elem_counts}
  {log}
  {color_scale}
  colorbar={{
    title: `MPtrj element counts by ${count_mode}`,
    title_style: `font-size: 1.3em;`,
  }}
/>
