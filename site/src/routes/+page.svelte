<script lang="ts">
  import MetricsTable from '$figs/metrics-table.svelte'
  import type { ModelData } from '$lib'
  import Readme from '$root/readme.md'
  import { onMount } from 'svelte'
  import all_stats from './models/model-stats.json'

  let best_model = Object.entries(all_stats).reduce((current, [model_name, stats]) => {
    if (!current?.F1 || stats.F1 > current.F1) {
      return { model_name, ...stats }
    }
    return current
  }, {}) as ModelData

  const metadata = import.meta.glob(`$root/models/**/metadata.yml`, {
    eager: true,
    import: `default`,
  }) as Record<string, ModelData | ModelData[]>

  onMount(async () => {
    if (best_model) {
      const md = metadata[`../models/${best_model.model_name.toLowerCase()}/metadata.yml`]
      best_model = { ...best_model, ...md }
    }
  })
</script>

<Readme>
  <div slot="best-report">
    {#if best_model}
      {@const { model_name, F1, R2, DAF, repo, doi } = best_model}
      We find <a href={repo}>{model_name}</a> (<a href={doi}>paper</a>) to achieve the
      highest F1 score of {F1}, R<sup>2</sup> of {R2}
      and a discovery acceleration factor (DAF) of {DAF}
      (meaning a ~{Number(DAF).toFixed(0)}x higher rate of stable structures compared to
      dummy selection in our already enriched search space).
    {/if}
  </div>
  <MetricsTable slot="metrics-table" />
</Readme>
