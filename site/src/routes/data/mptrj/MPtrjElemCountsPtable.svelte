<script lang="ts">
  import { PtableHeatmap } from '$lib'
  import figshare_urls from '$pkg/figshare/1.0.0.json'

  export let log = false // log color scale
  export let count_mode = `occurrence`

  const elem_counts = import.meta.glob(`../mp-trj-element-counts-by-*.json`, {
    eager: true,
    import: `default`,
  })

  $: mp_trj_elem_counts = elem_counts[`../mp-trj-element-counts-by-${count_mode}.json`]
</script>

<p>
  Element counts for
  <a href={figshare_urls.mptrj.article}>MPtrj training set</a> consisting of 1,580,395 structures
  which are frames of the DFT relaxations performed on all 154,719 MP materials.
</p>

<PtableHeatmap
  heatmap_values={mp_trj_elem_counts}
  {log}
  {...$$restProps}
  color_bar_props={{ label: `MPtrj element counts by ${count_mode}` }}
/>
