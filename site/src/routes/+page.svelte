<script lang="ts">
  import type { ModelData } from '$lib'
  import { CaptionedMetricsTable, model_is_compliant, MODEL_METADATA } from '$lib'
  import Readme from '$root/readme.md'
  import KappaNote from '$site/src/routes/kappa-note.md'
  import Icon from '@iconify/svelte'
  import { Toggle, Tooltip } from 'svelte-zoo'

  let show_non_compliant = false
  let show_energy_only = false

  $: best_model = MODEL_METADATA.reduce((best, md: ModelData) => {
    const best_F1 = best.metrics?.discovery?.full_test_set?.F1 ?? 0
    const md_F1 = md.metrics?.discovery?.full_test_set?.F1 ?? 0
    if ((!best_F1 || md_F1 > best_F1) && (show_non_compliant || model_is_compliant(md))) {
      return md
    }
    return best
  }, {} as ModelData)
</script>

<Readme>
  <span slot="model-count"
    >{MODEL_METADATA.filter((md) => show_non_compliant || model_is_compliant(md))
      .length}&ensp;</span
  >

  <div slot="best-report">
    {#if best_model}
      {@const { model_name, F1, R2, DAF, repo, paper } = best_model}
      {@const model_key = model_name.replaceAll(` `, `-`).toLowerCase()}
      <a href="/models/{model_key}">{model_name}</a> (<a href={paper}>paper</a>,
      <a href={repo}>code</a>) achieves the highest F1 score of {F1}, R<sup>2</sup> of {R2}
      and a discovery acceleration factor (DAF) of {DAF}
      (i.e. a ~{Number(DAF).toFixed(1)}x higher rate of stable structures compared to
      dummy discovery in the already enriched test set containing 16% stable materials).
    {/if}
  </div>

  <div slot="metrics-table" style="display: grid; gap: 1ex; place-items: center;">
    <KappaNote />
    <div style="display: flex; gap: 1em; align-items: center;">
      <Toggle bind:checked={show_non_compliant}
        >Show non-compliant models <Tooltip max_width="10em">
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
      <Toggle bind:checked={show_energy_only}
        >Show energy-only models <Tooltip max_width="10em">
          <span slot="tip">
            Models that only predict energy (E) perform worse<br /> and can't be evaluated
            on force-modeling tasks such as κ<sub>SRME</sub>
          </span>
          <Icon icon="octicon:info-16" inline style="padding: 0 3pt;" />
        </Tooltip>&ensp;</Toggle
      >
    </div>
    <CaptionedMetricsTable {show_non_compliant} {show_energy_only} />
  </div>
</Readme>
