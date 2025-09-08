<script lang="ts">
  import { MetricScatter, MetricsTable } from '$lib'
  import { ALL_METRICS, HYPERPARAMS, METADATA_COLS } from '$lib/labels'
  import KappaNote from './kappa-note.md'

  // Default column visibility
  let visible_cols: Record<string, boolean> = $state({
    // Hide other metrics
    ...Object.fromEntries(
      Object.values(ALL_METRICS).map((col) => [col.label, false]),
    ),
    // Show all metadata
    ...Object.fromEntries(
      Object.values(METADATA_COLS).map((col) => [col.label, true]),
    ),
    // Show phonon metrics
    [ALL_METRICS.κ_SRME.label]: true,
  })
</script>

<h1>MLFF Phonon Modeling Metrics</h1>

<section class="full-bleed">
  <MetricsTable col_filter={(col) => visible_cols[col.label] ?? true} sort_hint="" />
</section>

<KappaNote />

<h3>κ<sub>SRME</sub> vs Model Parameters</h3>
<p>
  κ<sub>SRME</sub> ranges from 0 to 2, the lower the better. This metric was introduced in
  <a href="https://arxiv.org/abs/2408.00755v4">arXiv:2408.00755v4</a>. This modeling task
  would not have been possible without the
  <a href="https://github.com/atztogo/phonondb">PhononDB</a>
  and the help of Atsushi Togo who kindly shared the
  <a
    href="https://github.com/atztogo/phonondb/blob/main/README.md#url-links-to-phono3py-finite-displacement-method-inputs-of-103-compounds-on-mdr-at-nims-pbe"
  >PBE reference data for the 103 MP structures</a> that form the test set for this task.
</p>

<MetricScatter
  x_prop={HYPERPARAMS.model_params}
  y_prop={ALL_METRICS.κ_SRME}
  style="height: 400px"
/>
