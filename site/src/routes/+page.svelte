<script lang="ts">
  import CumulativeMetricsPlot from '$root/figures/2022-12-26-cumulative-clf-metrics.svelte'
  import Readme from '$root/readme.md'
</script>

<Readme />

{#if typeof document !== `undefined`}
  <CumulativeMetricsPlot
    style="width: 90vw; transform: translateX(-45vw) translateX(min(25em, 45vw))"
  />
{/if}
