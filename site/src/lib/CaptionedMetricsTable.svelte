<script lang="ts">
  import { DiscoveryMetricsTable } from '$lib'
  import { pretty_num } from 'elementari'
  import { onMount } from 'svelte'

  export let show_non_compliant = false

  let n_wbm_stable_uniq_protos = 32_942
  let n_wbm_uniq_protos = 215_488

  onMount(() => {
    // convert model names into clickable links to /models/<model_key>
    for (const cell of document.querySelectorAll(`span[data-model-key]`)) {
      const model_key = cell.getAttribute(`data-model-key`)
      cell.innerHTML = `<a href="/models/${model_key}">${cell.innerHTML}</a>`
    }
  })
</script>

<figure {...$$props}>
  <DiscoveryMetricsTable {show_non_compliant} {...$$restProps} />
  <div class="downloads">
    Download table as
    {#each [`PDF`, `SVG`] as file_ext}
      {@const suffix = show_non_compliant ? `` : `-only-compliant`}
      <a href="/figs/metrics-table-uniq-protos{suffix}.{file_ext.toLowerCase()}" download>
        {file_ext}
      </a>
    {/each}
  </div>
  <figcaption>
    Training size is the number of materials used to train the model. For models trained
    on DFT relaxations, we show the number of distinct frames in parentheses. In cases
    where only the number of frames is known, we report the number of frames as the
    training set size. <code>(N=x)</code> in the Model Params column shows the number of
    estimators if an ensemble was used. DAF = Discovery Acceleration Factor measures how
    many more stable materials a model finds compared to random selection from the test
    set. The unique structure prototypes in the WBM test set have a
    <code>{pretty_num(n_wbm_stable_uniq_protos / n_wbm_uniq_protos, `.1%`)}</code> rate of
    stable crystals, meaning the max possible DAF is
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
  }
  figcaption {
    font-size: 0.9em;
    padding: 2pt 6pt;
    background-color: rgba(255, 255, 255, 0.07);
  }
  div.downloads {
    display: flex;
    gap: 1ex;
    justify-content: center;
    margin: 1ex 0;
  }
  div.downloads a {
    background-color: rgba(255, 255, 255, 0.1);
    padding: 0 6pt;
    border-radius: 4pt;
  }
</style>
