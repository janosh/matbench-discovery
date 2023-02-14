<script lang="ts">
  import type { ChemicalElement } from 'elementari'
  import { pretty_num } from 'elementari/labels'

  export let element: ChemicalElement
  export let elem_counts: number[] | Record<string, number>
  export let style: string | null = null

  $: count = elem_counts[element?.symbol] ?? elem_counts[element?.number] ?? null
  $: total = (
    Array.isArray(elem_counts) ? elem_counts : Object.values(elem_counts)
  ).reduce((a, b) => a + b, 0)
</script>

<strong {style}>
  {#if element?.name}
    {element?.name}: {pretty_num(count)}
    <!-- compute percent of total -->
    {#if count > 0}
      ({pretty_num((count / total) * 100)}%)
    {/if}
  {/if}
</strong>

<style>
  strong {
    display: block;
    text-align: center;
  }
</style>
