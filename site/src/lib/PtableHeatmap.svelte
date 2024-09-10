<script lang="ts">
  import { PtableInset } from '$lib'
  import { ColorBar, PeriodicTable, TableInset, type ChemicalElement } from 'elementari'
  import type { ComponentProps } from 'svelte'
  import { Toggle } from 'svelte-zoo'
  import type { Snapshot } from './$types'

  export let heatmap_values: Record<string, number>
  export let color_scale: string = `Viridis`
  export let active_element: ChemicalElement
  export let log = false // log color scale
  export let color_bar_props: ComponentProps<ColorBar> = {}

  export const snapshot: Snapshot = {
    capture: () => ({ color_scale, log }),
    restore: (values) => ({ color_scale, log } = values),
  }
</script>

<PeriodicTable
  {heatmap_values}
  {color_scale}
  {log}
  bind:active_element
  show_photo={false}
>
  <TableInset slot="inset">
    <label for="log">Log color scale<Toggle id="log" bind:checked={log} /></label>
    <PtableInset element={active_element} elem_counts={heatmap_values} />
    <ColorBar
      label="Count"
      label_side="top"
      {color_scale}
      tick_labels={5}
      range={[0, Math.max(...Object.values(heatmap_values))]}
      style="width: 85%; margin: 0 2em 2em;"
      {...color_bar_props}
    />
  </TableInset>
</PeriodicTable>

<style>
  label {
    display: flex;
    gap: 1ex;
    place-content: center;
    align-items: start;
    justify-items: center;
  }
</style>
