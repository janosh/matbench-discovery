<script lang="ts">
  import Readme from '$root/readme.md'
  import CumulativeMetricsPlot from '$site/figures/2022-12-26-cumulative-clf-metrics.svelte'
</script>

<Readme />

<br />

{#if typeof document !== `undefined`}
  <CumulativeMetricsPlot
    style="width: 90vw; transform: translateX(-45vw) translateX(min(25em, 45vw))"
  />
{/if}
