<script lang="ts">
  import { Icon } from '$lib'
  import type { IconName } from './icons'

  interface Props {
    icons: { id?: string; src?: string; name: string }[] | undefined
    [key: string]: unknown
  }
  let { icons = $bindable([]), ...rest }: Props = $props()
</script>

{#each icons ?? [] as logo (logo.id ?? logo.src)}
  {#if logo.id?.startsWith(`icon:`)}
    <span title={logo.name} {...rest}>
      <Icon icon={logo.id.replace(`icon:`, ``) as IconName} />
    </span>
  {:else if logo.src}
    <img
      src={logo.src}
      alt="{logo.name} logo"
      title={logo.name}
      {...rest}
      style={`margin: 0; height: 1em; ${rest.style ?? ``}`}
    />
  {/if}
{/each}

<style>
  :root[style*='color-scheme: light'] {
    --logo-brightness: 0.5;
  }
  span, img {
    filter: grayscale(100%) brightness(var(--logo-brightness, 1));
    height: 1em;
    width: auto;
    vertical-align: middle;
    margin: 0 !important;
  }
</style>
