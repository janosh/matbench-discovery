<script lang="ts">
  import type { OrgLogo } from '$lib'
  import { Icon } from 'svelte-widgets'

  // Suppress native titles when a parent supplies a richer tooltip.
  let { logo, show_title = true }: { logo: OrgLogo; show_title?: boolean } = $props()
</script>

{#if logo.icon}
  <span title={show_title ? logo.name : undefined} class="org-logo">
    <Icon icon={logo.icon} aria-label="{logo.name} logo" />
  </span>
{:else}
  <img
    src={logo.src}
    alt="{logo.name} logo"
    title={show_title ? logo.name : undefined}
    class="org-logo"
  />
{/if}

<style>
  .org-logo {
    filter: grayscale(100%) brightness(var(--logo-brightness, 1));
    height: 1em;
    width: auto;
    vertical-align: middle;
    margin: 0;
  }
  :global(:root[data-theme='light']) {
    --logo-brightness: 0.5;
  }
</style>
