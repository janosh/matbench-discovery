<script lang="ts">
  import type { DiscoverySet, ModelData } from '$lib/types'
  import { HeatmapTable } from 'matterviz'
  import type { Label, RowData } from 'matterviz'

  let {
    models,
    discovery_set,
  }: {
    models: ModelData[]
    discovery_set: DiscoverySet
  } = $props()

  const columns: Label[] = [
    { label: `Model`, sticky: true },
    {
      label: `Missing predictions`,
      description: `Fraction of the selected discovery set without a model prediction`,
      better: `lower`,
      format: `.2~%`,
    },
    {
      label: `False discoveries`,
      description: `Fraction of predicted-stable structures that are unstable according to DFT (1 − precision)`,
      better: `lower`,
      format: `.2~%`,
    },
    {
      label: `Missed stable`,
      description: `Fraction of DFT-stable structures predicted to be unstable (1 − recall)`,
      better: `lower`,
      format: `.2~%`,
    },
  ]

  const ratio = (numerator: number, denominator: number): number | null =>
    denominator > 0 ? numerator / denominator : null

  let table_data = $derived(
    models.flatMap((model): RowData[] => {
      const discovery = model.metrics?.discovery
      const metrics =
        discovery != null && typeof discovery === `object`
          ? discovery[discovery_set]
          : undefined
      if (!metrics) return []
      const classified = metrics.TP + metrics.FP + metrics.TN + metrics.FN
      return [
        {
          Model: `<a href="/models/${model.model_key}">${model.model_name}</a>`,
          // Missing predictions are filled as unstable before computing FN/TN, so
          // they are already included in the confusion-matrix total.
          'Missing predictions': ratio(metrics.missing_preds, classified),
          'False discoveries': ratio(metrics.FP, metrics.TP + metrics.FP),
          'Missed stable': ratio(metrics.FN, metrics.TP + metrics.FN),
        },
      ]
    }),
  )
</script>

<HeatmapTable
  data={table_data}
  {columns}
  initial_sort={{ column: `False discoveries`, direction: `asc` }}
  sort_hint=""
/>
