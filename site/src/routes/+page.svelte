<script lang="ts">
  import { dev } from '$app/environment'
  import type { ModelData, ModelStats } from '$lib'
  import { CaptionedMetricsTable } from '$lib'
  import all_stats from '$lib/model-stats.json'
  import Readme from '$root/readme.md'
  import { Toggle } from 'svelte-zoo'

  let show_proprietary = false

  const metadata = import.meta.glob(`$root/models/**/*.yml`, {
    eager: true,
    import: `default`,
  }) as Record<string, ModelData>

  $: best_model = Object.values(metadata).reduce((best, md: ModelData) => {
    const stats = all_stats[md.model_name] as ModelStats | undefined
    if (!stats && dev) {
      const avail_keys = Object.keys(all_stats)
      console.warn(
        `No metadata found for '${md.model_name}', available keys: ${avail_keys}`,
      )
      return best
    }

    const openness = md.open ?? `OSOD`
    if ((!best?.F1 || stats?.F1 > best?.F1) && (show_proprietary || openness == `OSOD`))
      return { ...stats, ...md } as ModelData

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
    <Toggle bind:checked={show_proprietary}>Show proprietary models ensp;</Toggle>
    <CaptionedMetricsTable bind:show_proprietary />
  </div>
</Readme>
