<script lang="ts">
  import { DISCOVERY_SET_LABELS, discovery_set_toggle_options } from '$lib/labels'
  import type { DiscoverySet } from '$lib/types'
  import { ButtonGroup, Icon } from 'svelte-widgets'
  import { Info } from 'svelte-widgets/icons'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    selected = $bindable(),
    ...rest
  }: { selected: DiscoverySet } & HTMLAttributes<HTMLDivElement> = $props()
</script>

<ButtonGroup
  bind:selected
  label="Discovery test set"
  options={discovery_set_toggle_options}
  tooltip_options={{ allow_html: true }}
  {...rest}
>
  {#snippet option_suffix({ option })}
    {@const link = DISCOVERY_SET_LABELS[option.value as DiscoverySet]?.link}
    {#if link}
      <a
        href={link}
        target="_blank"
        rel="noopener noreferrer"
        aria-label="About {option.label ?? option.value}"
        style="display: inline-flex; align-self: stretch; align-items: center; padding-right: 8px"
      >
        <Icon icon={Info} />
      </a>
    {/if}
  {/snippet}
</ButtonGroup>
