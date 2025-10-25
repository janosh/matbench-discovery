<script lang="ts">
  import { PtableInset } from '$lib'
  import {
    type ChemicalElement,
    ColorBar,
    PeriodicTable,
    TableInset,
  } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'
  import type { ComponentProps } from 'svelte'

  let {
    heatmap_values,
    color_scale = $bindable(`interpolateViridis`),
    active_element = $bindable(null),
    log = $bindable(false),
    colorbar = {},
    ...rest
  }: ComponentProps<typeof PeriodicTable> & {
    heatmap_values: Record<string, number>
    color_scale?: D3InterpolateName
    active_element?: ChemicalElement | null
    log?: boolean // log color scale
    colorbar?: ComponentProps<typeof ColorBar>
  } = $props()

  export const snapshot = {
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
  missing_color="rgba(255,255,255,0.3)"
  {...rest}
>
  {#snippet inset()}
    {@const style = `height: 1.5em; visibility: ${active_element ? `visible` : `hidden`};`}
    <TableInset>
      <label for="log">
        Log color scale<input id="log" type="checkbox" bind:checked={log} />
      </label>
      <PtableInset
        element={active_element}
        elem_counts={heatmap_values}
        {style}
      />
      <ColorBar
        title="Count"
        title_side="top"
        {color_scale}
        tick_labels={5}
        range={[0, Math.max(...Object.values(heatmap_values))]}
        style="width: 85%; margin: 0 2em 2em"
        {...colorbar}
      />
    </TableInset>
  {/snippet}
</PeriodicTable>

<style>
  label {
    display: flex;
    font-size: 1.1em;
    gap: 1ex;
    place-content: center;
    place-items: center;
  }
</style>
