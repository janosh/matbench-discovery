<script lang="ts">
  import { comparison } from '$lib/model-comparison.svelte'
  import { Icon } from 'svelte-widgets'
  import { Check, Scale } from 'svelte-widgets/icons'
  import type { HTMLButtonAttributes } from 'svelte/elements'

  let {
    model_key,
    model_name = model_key,
    open_dialog = false,
    ...rest
  }: HTMLButtonAttributes & {
    model_key: string
    model_name?: string
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
  aria-pressed={open_dialog ? undefined : selected}
  aria-label={open_dialog ? undefined : `${label} ${model_name}`}
  {title}
  {onclick}
  {...rest}
>
  <Icon icon={selected && !open_dialog ? Check : Scale} />
  {label}
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
