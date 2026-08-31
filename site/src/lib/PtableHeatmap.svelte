<script lang="ts">
  import { PtableInset } from '$lib'
  import { extent } from 'd3-array'
  import type { ChemicalElement, ElementSymbol } from 'matterviz'
  import { ColorBar, PeriodicTable, TableInset } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'
  import type { ComponentProps } from 'svelte'

  let {
    heatmap_values,
    color_scale = $bindable(`interpolateViridis`),
    active_element = $bindable(null),
    log = $bindable(false),
    colorbar = {},
    ...rest
  }: Omit<ComponentProps<typeof PeriodicTable>, `heatmap_values`> & {
    heatmap_values: Record<ElementSymbol, number>
    color_scale?: D3InterpolateName
    active_element?: ChemicalElement | null
    log?: boolean // Log color scale
    colorbar?: ComponentProps<typeof ColorBar>
  } = $props()

  // legend span mirrors the tile color ramp: PeriodicTable spans the finite data extent
  // and drops non-positive values in log mode (they have no log image)
  let colorbar_range = $derived.by((): [number, number] => {
    const finite_values = Object.values(heatmap_values).filter(
      (value) => Number.isFinite(value) && (!log || value > 0),
    )
    const [min, max] = extent(finite_values)
    return [min ?? 0, max ?? 0]
  })
</script>

<PeriodicTable
  {heatmap_values}
  {color_scale}
  {log}
  bind:active_element
  show_photo={false}
  missing={{ color: `rgba(255,255,255,0.3)` }}
  {...rest}
>
  {#snippet inset()}
    {@const style = `height: 1.5em; visibility: ${active_element ? `visible` : `hidden`};`}
    <TableInset>
      <label>
        Log color scale<input type="checkbox" bind:checked={log} />
      </label>
      {#if active_element}
        <PtableInset element={active_element} elem_counts={heatmap_values} {style} />
      {:else}
        <span {style}></span>
      {/if}
      <ColorBar
        title="Count"
        title_side="top"
        scale={color_scale}
        scale_type={log ? `log` : `linear`}
        tick_labels={5}
        range={colorbar_range}
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
