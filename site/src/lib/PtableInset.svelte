<script lang="ts">
  import type { ChemicalElement } from 'elementari'
  import { pretty_num } from 'elementari'

  export let element: ChemicalElement
  export let elem_counts: number[] | Record<string, number>
  export let style: string | null = null
  export let show_percent: boolean = true
  export let unit: string = ``

  $: value = elem_counts[element?.symbol] ?? elem_counts[element?.number - 1] ?? null
  $: total = (
    Array.isArray(elem_counts) ? elem_counts : Object.values(elem_counts)
  ).reduce((a, b) => a + b, 0)
</script>

<strong {style}>
  {#if element?.name}
    {element?.name}: {pretty_num(value)}
    {@html unit}
    <!-- compute percent of total -->
    {#if show_percent && total > 0}
      ({pretty_num((100 * value) / total)}%)
    {/if}
  {/if}
</strong>

<style>
  strong {
    display: block;
    text-align: center;
    min-height: 18pt;
  }
</style>
