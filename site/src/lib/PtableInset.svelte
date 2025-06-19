<script lang="ts">
  import type { ChemicalElement, ElementSymbol } from 'matterviz'
  import { format_num } from 'matterviz'

  interface Props {
    element: ChemicalElement
    elem_counts: number[] | Record<ElementSymbol, number>
    show_percent?: boolean
    unit?: string
    [key: string]: unknown
  }
  let { element, elem_counts, show_percent = true, unit = ``, ...rest }: Props = $props()

  let value = $derived(
    Array.isArray(elem_counts)
      ? (elem_counts[element?.number - 1] ?? null)
      : (elem_counts[element?.symbol] ?? null),
  )
  let total = $derived(
    (Array.isArray(elem_counts) ? elem_counts : Object.values(elem_counts)).reduce(
      (total, count) => total + count,
      0,
    ),
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
