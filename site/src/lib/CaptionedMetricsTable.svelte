<script lang="ts">
  import MetricsTable from '$figs/metrics-table.svelte'
</script>

<figure {...$$props}>
  <MetricsTable />
  <figcaption>
    Training size is the number of structures used to train the model. For models trained
    on DFT relaxations, we show the number of distinct materials in parentheses where we
    count multiple frames from the trajectory as one material.
  </figcaption>
</figure>

<style>
  figure {
    margin: 0;
  }
  figcaption {
    font-size: 0.9em;
    text-align: center;
    margin-top: 1ex;
  }
</style>
