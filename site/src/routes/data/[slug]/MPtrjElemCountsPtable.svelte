<script lang="ts">
  import { data_files, PtableHeatmap } from '$lib'
  import type { ComponentProps } from 'svelte'

  interface Props extends Partial<ComponentProps<typeof PtableHeatmap>> {
    count_mode?: string
  }
  let { count_mode = `occurrence`, ...rest }: Props = $props()

  const elem_counts = import.meta.glob(`../mp-trj-element-counts-by-*.json`, {
    eager: true,
    import: `default`,
  }) as Record<string, Record<string, number>>

  let mp_trj_elem_counts =
    elem_counts[`../mp-trj-element-counts-by-${count_mode}.json`]

  if (!(`mp_trj_json_gz` in data_files)) {
    throw `mp_trj_json_gz not found in data-files.yml`
  }
</script>

<p>
  Element counts for
  <a href={data_files.mp_trj_json_gz.figshare}>MPtrj training set</a> consisting of
  1,580,395 structures which are frames of the DFT relaxations performed on all 154,719 MP
  materials.
</p>

<PtableHeatmap
  {...rest}
  heatmap_values={mp_trj_elem_counts}
  colorbar={{
    title: `MPtrj element counts by ${count_mode}`,
    title_style: `font-size: 1.3em;`,
  }}
/>
