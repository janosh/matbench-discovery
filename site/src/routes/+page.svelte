<script lang="ts">
  import CumulativeMetricsPlot from '$root/figures/2022-12-25-precision-recall-curves.svelte'
  import Readme from '$root/readme.md'
  import { onMount } from 'svelte'

  onMount(() => {
    // make links to site/src/routes/ files deployment-compatible
    const site_links: NodeListOf<HTMLAnchorElement> = document.querySelectorAll(
      `a[href^='site/src/routes']`
    )
    for (const link of site_links) {
      link.href = link.href.replace(`/site/src/routes`, ``).split('/+page')[0]
    }
  })
</script>

<Readme />

{#if typeof document !== `undefined`}
  <CumulativeMetricsPlot
    style="width: 90vw; transform: translateX(-45vw) translateX(min(25em, 45vw))"
  />
{/if}
