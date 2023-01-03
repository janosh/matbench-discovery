<script lang="ts">
  import FormEnergyHist from '$root/data/wbm/2022-12-07-hist-e-form-per-atom.svelte'
  import DataReadme from '$root/data/wbm/readme.md'
  // import { PeriodicTable, Toggle } from 'periodic-tables'
  import PeriodicTable from 'periodic-tables'
  import Toggle from 'periodic-tables/Toggle.svelte'
  import elem_counts from './2022-12-30-wbm-element-counts.json'

  const heatmap_values = Object.values(elem_counts)
  const color_map = {
    200: `blue`,
    35_000: `green`,
    80_000: `yellow`,
    150_000: `red`,
  }

  let log = true
</script>

<DataReadme>
  <svelte:fragment slot="hist-e-form-per-atom">
    {#if typeof document !== `undefined`}
      <FormEnergyHist />
    {/if}
  </svelte:fragment>
  <svelte:fragment slot="wbm-elements-log">
    <span>Log color scale? <Toggle bind:checked={log} /></span>
    <PeriodicTable {heatmap_values} {color_map} {log} />
  </svelte:fragment>
</DataReadme>

<style>
  span {
    display: flex;
    gap: 1ex;
    position: absolute;
    left: 50%;
    transform: translateX(-50%);
    z-index: 1;
  }
</style>
