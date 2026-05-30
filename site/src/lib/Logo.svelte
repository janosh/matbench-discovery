<script lang="ts">
  import type { OrgLogo } from '$lib/labels'
  import { Icon } from 'matterviz'

  // show_title controls the native title attribute. Disable it when the logo is
  // rendered inside another element that already provides a richer tooltip to
  // avoid duplicate/competing tooltips.
  let { logo, show_title = true }: { logo: OrgLogo; show_title?: boolean } = $props()
</script>

{#if logo.validated_icon}
  <span title={show_title ? logo.name : undefined} class="org-logo">
    <Icon icon={logo.validated_icon} />
  </span>
{:else if logo.src}
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
