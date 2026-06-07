<script lang="ts">
  import type { ChemicalElement, ElementSymbol } from 'matterviz'
  import { format_num } from 'matterviz'
  import type { HTMLAttributes } from 'svelte/elements'
  import { sum } from 'd3-array'

  let { element, elem_counts, show_percent = true, unit = ``, ...rest }: {
    element: ChemicalElement
    elem_counts: number[] | Record<ElementSymbol, number>
    show_percent?: boolean
    unit?: string
  } & HTMLAttributes<HTMLElementTagNameMap[`strong`]> = $props()

  let value = $derived(
    Array.isArray(elem_counts)
      ? (elem_counts[element?.number - 1] ?? null)
      : (elem_counts[element?.symbol] ?? null),
  )
  let total = $derived(
    sum(Array.isArray(elem_counts) ? elem_counts : Object.values(elem_counts)),
  )
</script>

<strong {...rest}>
  {#if element?.name}
    {element?.name}: {format_num(value)}
    {@html unit}
    <!-- compute percent of total -->
    {#if show_percent && total > 0}
      ({format_num((100 * value) / total)}%)
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
