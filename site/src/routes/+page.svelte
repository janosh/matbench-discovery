<script lang="ts">
  import { browser } from '$app/environment'
  import CumulativeMetricsPlot from '$figs/2023-01-31-cumulative-clf-metrics.svelte'
  import Readme from '$root/readme.md'
</script>

<Readme />

<br />

{#if browser}
  <CumulativeMetricsPlot
    style="width: 90vw; transform: translateX(-45vw) translateX(min(25em, 45vw))"
  />
{/if}
