<script lang="ts">
  import { PtableInset } from '$lib'
  import { ColorBar, PeriodicTable, TableInset, type ChemicalElement } from 'elementari'
  import type { ComponentProps } from 'svelte'
  import { Toggle } from 'svelte-zoo'
  import type { Snapshot } from './$types'

  interface Props {
    heatmap_values: Record<string, number>
    color_scale?: string
    active_element: ChemicalElement | null
    log?: boolean // log color scale
    color_bar_props?: ComponentProps<typeof ColorBar>
  }
  let {
    heatmap_values,
    color_scale = $bindable(`Viridis`),
    active_element = $bindable(null),
    log = $bindable(false),
    color_bar_props = {},
  }: Props = $props()

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
  {#snippet inset()}
    <TableInset>
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
  {/snippet}
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
