<script lang="ts">
  import Plot from '$root/figures/2022-12-25-precision-recall-curves.svelte'
  import Readme from '$root/readme.md'
</script>

<Readme />

{#if typeof document !== `undefined`}
  <Plot style="margin: 0 -5vw;" />
{/if}
