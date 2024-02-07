<script lang="ts">
  import MetricsTable from '$figs/metrics-table-uniq-protos.svelte'

  export let show_proprietary = false
</script>

<figure {...$$props} class:hide-prop={!show_proprietary}>
  <MetricsTable />
  <figcaption>
    Training size is the number of structures used to train the model. For models trained
    on DFT relaxations, we show the number of distinct materials in parentheses where we
    count multiple frames from the same trajectory as one material.
  </figcaption>
</figure>

<style>
  figure {
    margin: 0;
    display: grid;
    gap: 1ex;
    overflow: scroll;
    max-width: 90vw; /* enable horizontal scrolling on smaller screens */
  }
  figcaption {
    font-size: 0.9em;
    padding: 2pt 6pt;
    background-color: rgba(255, 255, 255, 0.07);
  }
  /* hide rows (<tr>) where any cell has a class of proprietary */
  figure.hide-prop :global(tr:has(.proprietary)) {
    display: none;
  }
</style>
