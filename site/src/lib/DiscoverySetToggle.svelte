<script lang="ts">
  import { DISCOVERY_SET_LABELS, discovery_set_toggle_options } from '$lib/labels'
  import type { DiscoverySet } from '$lib/types'
  import { ButtonGroup, Icon } from 'svelte-widgets'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    selected = $bindable(),
    ...rest
  }: { selected: DiscoverySet } & HTMLAttributes<HTMLDivElement> = $props()
</script>

<ButtonGroup
  bind:selected
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
        style="display: inline-flex; align-items: center; padding-right: 8px; font-size: 0.83em"
      >
        <Icon icon="Info" style="transform: translateY(-1px)" />
      </a>
    {/if}
  {/snippet}
</ButtonGroup>
