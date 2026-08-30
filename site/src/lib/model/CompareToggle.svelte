<script lang="ts">
  import { comparison } from '$lib/model-comparison.svelte'
  import { Icon } from 'svelte-widgets'
  import { Check, Scale } from 'svelte-widgets/icons'
  import type { HTMLButtonAttributes } from 'svelte/elements'

  let { model_key, ...rest }: HTMLButtonAttributes & { model_key: string } = $props()
  let selected = $derived(comparison.keys.has(model_key))
</script>

<button
  class:selected
  aria-pressed={selected}
  title="{selected ? `Remove from` : `Add to`} model comparison"
  onclick={() => comparison.toggle(model_key)}
  {...rest}
>
  <Icon icon={selected ? Check : Scale} />
  {selected ? `Comparing` : `Compare`}
</button>

<style>
  button {
    display: inline-flex;
    align-items: center;
    gap: 5px;
    padding: 0 5pt;
    background: var(--chip-bg);
    border-radius: 5px;
    font: inherit;
    color: inherit;
  }
  button.selected {
    color: var(--link-color);
    box-shadow: inset 0 0 0 1px var(--link-color);
  }
</style>
