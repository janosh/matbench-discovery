<script lang="ts">
  import MetricsTable from '$figs/metrics-table-uniq-protos.svelte'
  import { pretty_num } from 'elementari'

  export let show_non_compliant = false

  let n_wbm_stable_uniq_protos = 32_942
  let n_wbm_uniq_protos = 215_488
</script>

<figure {...$$props} class:nide-non-compliant={!show_non_compliant}>
  <MetricsTable />
  <figcaption>
    Training size is the number of materials used to train the model. For models trained
    on DFT relaxations, we show the number of distinct frames in parentheses. In cases
    where only the number of frames is known, we report the number of frames as the
    training set size. <code>(N=x)</code> in the Model Params column shows the number of
    estimators if an ensemble was used. DAF = Discovery Acceleration Factor measures how
    many more stable materials a model finds compared to random selection from the test
    set. The WBM test set has a 16.7% rate of stable crystals, meaning the max possible
    DAF is
    <code
      >({pretty_num(n_wbm_stable_uniq_protos)} / {pretty_num(n_wbm_uniq_protos)})^−1 ≈
      {pretty_num(n_wbm_uniq_protos / n_wbm_stable_uniq_protos)}</code
    >.
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
  /* hide rows (<tr>) where any cell has a class of non-compliant */
  figure.nide-non-compliant :global(tr:has(.non-compliant)) {
    display: none;
  }
</style>
