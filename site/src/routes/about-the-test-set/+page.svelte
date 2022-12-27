<script lang="ts">
  import Plot from '$root/data/wbm/2022-12-07-hist-e-form-per-atom.svelte'
  import DataReadme from '$root/data/wbm/readme.md'
</script>

<DataReadme>
  <svelte:fragment slot="hist-e-form-per-atom">
    {#if typeof document !== `undefined`}
      <Plot />
    {/if}
  </svelte:fragment>
</DataReadme>
