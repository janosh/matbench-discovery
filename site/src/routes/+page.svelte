<script lang="ts">
  import CumulativeMetricsPlot from '$figs/2023-01-19-cumulative-clf-metrics.svelte'
  import Readme from '$root/readme.md'
</script>

<Readme />

<br />

{#if typeof document !== `undefined`}
  <CumulativeMetricsPlot
    style="width: 90vw; transform: translateX(-45vw) translateX(min(25em, 45vw))"
  />
{/if}
