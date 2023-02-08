<script lang="ts">
  import { browser } from '$app/environment'
  import CumulativeMetricsPlot from '$figs/cumulative-clf-metrics.svelte'
  import MetricsTable from '$figs/metrics-table.svelte'
  import Readme from '$root/readme.md'
</script>

<Readme>
  <div slot="metrics-table" style="container-type: inline-size;">
    <MetricsTable style="font-size: 1.65cqw;" />
  </div>
</Readme>

<br />

{#if browser}
  <CumulativeMetricsPlot
    style="width: 90vw; transform: translateX(-45vw) translateX(min(25em, 45vw))"
  />
{/if}
