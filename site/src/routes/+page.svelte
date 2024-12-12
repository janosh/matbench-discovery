<script lang="ts">
  import type { ModelData } from '$lib'
  import { CaptionedMetricsTable, model_is_compliant, MODEL_METADATA } from '$lib'
  import Readme from '$root/readme.md'
  import KappaNote from '$site/src/routes/kappa-note.md'
  import Icon from '@iconify/svelte'
  import { Toggle, Tooltip } from 'svelte-zoo'

  let show_non_compliant: boolean = false
  let show_energy_only: boolean = false

  // Default column visibility
  let visible_cols: Record<string, boolean> = {
    Model: true,
    F1: true,
    DAF: true,
    Prec: true,
    Acc: true,
    TPR: false,
    TNR: false,
    MAE: true,
    RMSE: false,
    'R<sup>2</sup>': true,
    'κ<sub>SRME</sub>': true,
    'Training Set': true,
    Params: true,
    Targets: true,
    'Date Added': true,
  }

  $: best_model = MODEL_METADATA.reduce((best, md: ModelData) => {
    const best_F1 = best.metrics?.discovery?.full_test_set?.F1 ?? 0
    const md_F1 = md.metrics?.discovery?.full_test_set?.F1 ?? 0
    if ((!best_F1 || md_F1 > best_F1) && (show_non_compliant || model_is_compliant(md))) {
      return md
    }
    return best
  }, {} as ModelData)

  // Get array of hidden columns
  $: hide_cols = Object.entries(visible_cols)
    .filter(([_, visible]) => !visible)
    .map(([col]) => col)

  let column_panel_open: boolean = false
</script>

<svelte:body
  on:click={(event) => {
    if (!event.target?.closest(`.column-toggles`)) {
      column_panel_open = false
    }
  }}
/>

<Readme>
  <CaptionedMetricsTable
    {show_non_compliant}
    {show_energy_only}
    {hide_cols}
    style="margin-top: 4em;"
    slot="metrics-table"
  />
  <div slot="table-controls">
    <Toggle bind:checked={show_non_compliant} style="gap: 3pt;">
      Show non-compliant models <Tooltip max_width="20em">
        <span slot="tip">
          Models can be non-compliant for multiple reasons<br />
          - closed source (model implementation and/or train/test code)<br />
          - closed weights<br />
          - trained on more than the permissible training set (<a
            href="https://docs.materialsproject.org/changes/database-versions#v2022.10.28"
            >MP v2022.10.28 release</a
          >)<br />
          We still show these models behind a toggle as we expect them<br /> to nonetheless
          provide helpful signals for developing future models.
        </span>
        <Icon icon="octicon:info-16" inline style="padding: 0 3pt;" />
      </Tooltip>&ensp;</Toggle
    >
    <Toggle bind:checked={show_energy_only} style="gap: 3pt;">
      Show energy-only models <Tooltip max_width="20em">
        <span slot="tip">
          Models that only predict energy (E) perform worse<br /> and can't be evaluated
          on force-modeling tasks such as κ<sub>SRME</sub>
        </span>
        <Icon icon="octicon:info-16" inline style="padding: 0 3pt;" />
      </Tooltip>&ensp;</Toggle
    >

    <details class="column-toggles" bind:open={column_panel_open}>
      <summary>
        Columns <Icon icon="octicon:columns-16" inline />
      </summary>
      <div class="column-menu">
        {#each Object.keys(visible_cols) as col}
          <label>
            <input type="checkbox" bind:checked={visible_cols[col]} />
            {@html col}
          </label>
        {/each}
      </div>
    </details>
  </div>

  <span slot="model-count">
    {MODEL_METADATA.filter((md) => show_non_compliant || model_is_compliant(md)).length}
  </span>

  <div slot="best-report">
    {#if best_model}
      {@const { model_name, model_key, F1, R2, DAF, repo, paper } = best_model}
      <a href="/models/{model_key}">{model_name}</a> (<a href={paper}>paper</a>,
      <a href={repo}>code</a>) achieves the highest F1 score of {F1}, R<sup>2</sup> of {R2}
      and a discovery acceleration factor (DAF) of {DAF}
      (i.e. a ~{Number(DAF).toFixed(1)}x higher rate of stable structures compared to
      dummy discovery in the already enriched test set containing 16% stable materials).
    {/if}
  </div>
</Readme>
<KappaNote />

<style>
  div[slot='table-controls'] {
    display: flex;
    flex-wrap: wrap;
    gap: 1em;
    place-content: center;
    margin: 2em auto;
    max-width: 500px;
  }
  .column-toggles {
    position: relative;
  }
  .column-toggles summary {
    background: rgba(255, 255, 255, 0.1);
    padding: 2pt 6pt;
    border-radius: 4pt;
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 4px;
  }
  .column-toggles summary:hover {
    background: rgba(255, 255, 255, 0.15);
  }
  .column-toggles summary::-webkit-details-marker {
    display: none;
  }
  .column-menu {
    position: absolute;
    right: 0;
    top: calc(100% + 4pt);
    background: #1c1c1c;
    border: 1px solid rgba(255, 255, 255, 0.1);
    border-radius: 4pt;
    padding: 0 4pt;
    min-width: 150px;
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(120px, 1fr));
  }
  .column-menu label {
    display: inline-block;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }
  .column-menu label:hover {
    background: rgba(255, 255, 255, 0.1);
  }
  .column-menu :global(:is(sub, sup)) {
    transform: translate(-3pt, 6pt);
    font-size: 0.7em;
  }
</style>
