<script lang="ts">
  import { browser } from '$app/environment'
  import CumulativeMetricsPlot from '$figs/cumulative-clf-metrics.svelte'
  import MetricsTable from '$figs/metrics-table.svelte'
  import Readme from '$root/readme.md'
</script>

<Readme>
  <MetricsTable slot="metrics-table" />
</Readme>

<br />

{#if browser}
  <CumulativeMetricsPlot
    style="width: 90vw; transform: translateX(-45vw) translateX(min(25em, 45vw))"
  />
{/if}
