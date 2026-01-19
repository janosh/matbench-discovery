<script lang="ts">
  import type { OrgLogo } from '$lib/labels'
  import { Icon } from 'matterviz'

  let { logo }: { logo: OrgLogo } = $props()
</script>

{#if logo.validated_icon}
  <span title={logo.name} class="org-logo">
    <Icon icon={logo.validated_icon} />
  </span>
{:else if logo.src}
  <img
    src={logo.src}
    alt="{logo.name} logo"
    title={logo.name}
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
