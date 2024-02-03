<script lang="ts">
  import { dev } from '$app/environment'
  import { CaptionedMetricsTable, type ModelData } from '$lib'
  import Readme from '$root/readme.md'
  import { Toggle } from 'svelte-zoo'
  import all_stats from './models/model-stats.json'

  let show_proprietary = false

  const metadata = import.meta.glob(`$root/models/**/metadata.yml`, {
    eager: true,
    import: `default`,
  }) as Record<string, ModelData | ModelData[]>

  $: best_model = Object.entries(all_stats).reduce((best, [model_name, stats]) => {
    let md = metadata[`../models/${model_name.toLowerCase()}/metadata.yml`]
    if (!md && dev) {
      console.warn(`No metadata found for ${model_name}`)
    }
    if (Array.isArray(md)) md = md[0]

    const model_data = { ...stats, ...md } as ModelData

    const openness = model_data.open ?? `OSOD`
    if (
      (!best?.F1 || model_data.F1 > best.F1) &&
      (show_proprietary || openness == `OSOD`)
    )
      return model_data

    return best
  }, {} as ModelData)
</script>

<Readme>
  <span slot="model-count">{Object.keys(all_stats).length}&ensp;</span>

  <div slot="best-report">
    {#if best_model}
      {@const { model_name, F1, R2, DAF, repo, doi } = best_model}
      We find <a href={repo}>{model_name}</a> (<a href={doi}>paper</a>) to achieve the
      highest F1 score of {F1}, R<sup>2</sup> of {R2}
      and a discovery acceleration factor (DAF) of {DAF}
      (i.e. a ~{Number(DAF).toFixed(1)}x higher rate of stable structures compared to
      dummy discovery in the already enriched test set containing 16% stable materials).
    {/if}
  </div>

  <div slot="metrics-table" style="display: grid; gap: 1ex; place-items: center;">
    <Toggle bind:checked={show_proprietary}>Show proprietary models</Toggle>
    <CaptionedMetricsTable bind:show_proprietary />
  </div>
</Readme>
