<script lang="ts">
  import { comparison } from '$lib/model-comparison.svelte'
  import { Icon } from 'svelte-widgets'
  import { Check, Scale } from 'svelte-widgets/icons'
  import type { HTMLButtonAttributes } from 'svelte/elements'

  let {
    model_key,
    model_name = model_key,
    compact = false,
    open_dialog = false,
    ...rest
  }: HTMLButtonAttributes & {
    model_key: string
    model_name?: string
    // icon-only variant for table rows: hidden until the row is hovered, unless selected
    compact?: boolean
    // model pages: always add (never remove) the model and open the comparison dialog,
    // so the click lands the user in the picker instead of silently changing a set
    open_dialog?: boolean
  } = $props()
  let selected = $derived(comparison.keys.has(model_key))
  let n_others = $derived(comparison.keys.size - (selected ? 1 : 0))
  let label = $derived.by(() => {
    if (!open_dialog) return selected ? `Comparing` : `Compare`
    if (selected && n_others > 0)
      return `Comparing with ${n_others} other${n_others === 1 ? `` : `s`}`
    return `Compare with…`
  })
  let title = $derived(
    open_dialog
      ? `Open the side-by-side comparison with this model`
      : `${selected ? `Remove from` : `Add to`} model comparison`,
  )
  const onclick = () =>
    open_dialog ? comparison.open_with(model_key) : comparison.toggle(model_key)
</script>

<button
  class:selected
  class:compact
  aria-pressed={open_dialog ? undefined : selected}
  aria-label={compact ? `${label} ${model_name}` : undefined}
  {title}
  {onclick}
  {...rest}
>
  <Icon icon={selected && !open_dialog ? Check : Scale} />
  {#if !compact}{label}{/if}
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
  /* row action: takes up its space always (no layout shift) but only shows on row hover,
     keyboard focus, when selected, or on touch devices which have no hover */
  button.compact {
    margin-left: 4pt;
    padding: 0 2pt;
    background: none;
    vertical-align: middle;
    opacity: 0;
    transition: opacity 0.15s;
  }
  :global(tr:is(:hover, :focus-within)) button.compact,
  button.compact:is(.selected, :focus-visible) {
    opacity: 1;
  }
  @media (hover: none) {
    button.compact {
      opacity: 1;
    }
  }
</style>
