<script lang="ts">
  import FormEnergyHist from '$root/data/wbm/2022-12-07-hist-e-form-per-atom.svelte'
  import DataReadme from '$root/data/wbm/readme.md'
  import type { ChemicalElement } from 'sveriodic-table'
  import { PeriodicTable, TableInset, Toggle } from 'sveriodic-table'
  import { pretty_num } from 'sveriodic-table/labels'
  import elem_counts from './2022-12-30-wbm-element-counts.json'

  let log_color_scale = false
  const heatmap_values: number[] = Object.values(elem_counts)
  const color_map = {
    200: `blue`,
    35_000: `green`,
    80_000: `yellow`,
    150_000: `red`,
  }
  let active_element: ChemicalElement
</script>

<DataReadme>
  <svelte:fragment slot="hist-e-form-per-atom">
    {#if typeof document !== `undefined`}
      <FormEnergyHist />
    {/if}
  </svelte:fragment>
  <svelte:fragment slot="wbm-elements-log">
    <span>Log color scale <Toggle bind:checked={log_color_scale} /></span>
    <PeriodicTable {heatmap_values} {color_map} log={log_color_scale} bind:active_element>
      <TableInset slot="inset" grid_row="3">
        {#if active_element?.name}
          <strong>
            {active_element?.name}: {pretty_num(elem_counts[active_element?.symbol])}
            <!-- compute percent of total -->
            {#if elem_counts[active_element?.symbol] > 0}
              {@const total = heatmap_values.reduce((a, b) => a + b, 0)}
              ({pretty_num((elem_counts[active_element?.symbol] / total) * 100)}%)
            {/if}
          </strong>
        {/if}
      </TableInset>
    </PeriodicTable>
  </svelte:fragment>
</DataReadme>

<style>
  span {
    display: flex;
    gap: 1ex;
    position: absolute;
    left: 50%;
    transform: translateX(-50%) translateY(100%);
    z-index: 1;
  }
  strong {
    text-align: center;
  }
</style>
