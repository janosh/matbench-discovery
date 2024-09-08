<script lang="ts">
  import { dev } from '$app/environment'
  import type { ModelData, ModelStats } from '$lib'
  import { CaptionedMetricsTable } from '$lib'
  import all_stats from '$lib/model-stats-uniq-protos.json'
  import Readme from '$root/readme.md'
  import Icon from '@iconify/svelte'
  import { Toggle, Tooltip } from 'svelte-zoo'

  let show_non_compliant = false

  const metadata = import.meta.glob(`$root/models/**/*.yml`, {
    eager: true,
    import: `default`,
  }) as Record<string, ModelData>

  $: best_model = Object.values(metadata).reduce((best, md: ModelData) => {
    const stats = all_stats[md.model_name] as ModelStats | undefined
    if (!stats && dev) {
      const avail_keys = Object.keys(all_stats).join(`, `)
      console.warn(`No metadata for '${md.model_name}', available keys: ${avail_keys}`)
      return best
    }

    const openness = md.openness ?? `OSOD`
    if ((!best?.F1 || stats?.F1 > best?.F1) && (show_non_compliant || openness == `OSOD`))
      return { ...stats, ...md } as ModelData

    return best
  }, {} as ModelData)
</script>

<Readme>
  <span slot="model-count">{Object.keys(all_stats).length}&ensp;</span>

  <div slot="best-report">
    {#if best_model}
      {@const { model_name, F1, R2, DAF, repo, paper, model_type } = best_model}
      {@const model_key = model_name.replaceAll(` `, `-`).toLowerCase()}
      <a href="/models/{model_key}">{model_name}</a> (<a href={paper}>paper</a>,
      <a href={repo}>code</a>) achieves the highest F1 score of {F1}, R<sup>2</sup> of {R2}
      and a discovery acceleration factor (DAF) of {DAF}
      (i.e. a ~{Number(DAF).toFixed(1)}x higher rate of stable structures compared to
      dummy discovery in the already enriched test set containing 16% stable materials).
    {/if}
  </div>

  <div slot="metrics-table" style="display: grid; gap: 1ex; place-items: center;">
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
    <CaptionedMetricsTable bind:show_non_compliant />
  </div>
</Readme>
