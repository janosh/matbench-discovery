<script lang="ts">
  import FormEnergyHist from '$root/data/wbm/2022-12-07-hist-e-form-per-atom.svelte'
  import DataReadme from '$root/data/wbm/readme.md'
  import type { ChemicalElement } from 'sveriodic-table'
  import { PeriodicTable, TableInset, Toggle } from 'sveriodic-table'
  import { pretty_num } from 'sveriodic-table/labels'
  import mp_elem_counts from './2023-01-08-mp-element-counts.json'
  import wbm_elem_counts from './2023-01-08-wbm-element-counts.json'

  let log = false // log color scale
  const wbm_heat_vals: number[] = Object.values(wbm_elem_counts)
  const mp_heat_vals: number[] = Object.values(mp_elem_counts)
  const color_map = {
    200: `blue`,
    35_000: `green`,
    80_000: `yellow`,
    150_000: `red`,
  }
  let active_mp_elem: ChemicalElement
  let active_wbm_elem: ChemicalElement
</script>

<DataReadme>
  <svelte:fragment slot="hist-e-form-per-atom">
    {#if typeof document !== `undefined`}
      <FormEnergyHist />
    {/if}
  </svelte:fragment>
  <svelte:fragment slot="wbm-elements-heatmap">
    <span>Log color scale <Toggle bind:checked={log} /></span>
    <PeriodicTable heatmap_values={wbm_heat_vals} {color_map} {log} bind:active_element={active_wbm_elem}>
      <TableInset slot="inset" grid_row="3">
        {#if active_wbm_elem?.name}
          <strong>
            {active_wbm_elem?.name}: {pretty_num(wbm_elem_counts[active_wbm_elem?.symbol])}
            <!-- compute percent of total -->
            {#if wbm_elem_counts[active_wbm_elem?.symbol] > 0}
              {@const total = wbm_heat_vals.reduce((a, b) => a + b, 0)}
              ({pretty_num((wbm_elem_counts[active_wbm_elem?.symbol] / total) * 100)}%)
            {/if}
          </strong>
        {/if}
      </TableInset>
    </PeriodicTable>
  </svelte:fragment>
  <svelte:fragment slot="mp-elements-heatmap">
    <PeriodicTable heatmap_values={mp_heat_vals} {color_map} {log} bind:active_element={active_mp_elem}>
      <TableInset slot="inset" grid_row="3">
        {#if active_mp_elem?.name}
          <strong>
            {active_mp_elem?.name}: {pretty_num(wbm_elem_counts[active_mp_elem?.symbol])}
            <!-- compute percent of total -->
            {#if wbm_elem_counts[active_mp_elem?.symbol] > 0}
              {@const total = wbm_heat_vals.reduce((a, b) => a + b, 0)}
              ({pretty_num((wbm_elem_counts[active_mp_elem?.symbol] / total) * 100)}%)
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
