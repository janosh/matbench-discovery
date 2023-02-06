<script lang="ts">
  import { browser } from '$app/environment'
  import Readme from '$root/readme.md'
  import CumulativeMetricsPlot from '$site/src/figs/cumulative-clf-metrics.svelte'
</script>

<Readme />

<br />

{#if browser}
  <CumulativeMetricsPlot
    style="width: 90vw; transform: translateX(-45vw) translateX(min(25em, 45vw))"
  />
{/if}
