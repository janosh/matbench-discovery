<script lang="ts">
  import box_data from '$figs/box-hull-dist-errors.jsonl'
  import cum_pr from '$figs/cumulative-precision-recall.jsonl'
  import hist_clf from '$figs/hist-clf-pred-hull-dist.jsonl'
  import roc from '$figs/roc-models.jsonl'
  import rolling_mae from '$figs/rolling-mae-vs-hull-dist.jsonl'
  import { dashed, labeled_vline, model_mae, order_models, wide_legend } from '$lib/fig-helpers'
  import { model_is_compliant, MODELS } from '$lib/models.svelte'
  import { BarPlot, BoxPlot, PlotLegend, ScatterPlot } from 'matterviz/plot'
  import type { DataSeries, FillRegion, LegendItem } from 'matterviz/plot'

  let { show_non_compliant = true }: { show_non_compliant?: boolean } = $props()

  // payload model entries carry `key` (= MODELS model_key) for this compliance join
  const compliant_keys = new Set(
    MODELS.filter(model_is_compliant).map((model) => model.model_key),
  )
  const shown = <T extends { key: string }>(models: T[]): T[] =>
    show_non_compliant ? models : models.filter((mdl) => compliant_keys.has(mdl.key))

  // payload models arrive pre-styled (stable MODELS colors + discovery-F1-desc leaderboard
  // order) from the figure_payload plugin. cumulative P/R keeps that default order; the rest
  // re-rank for their figure: box + rolling by discovery MAE (best first, rolling shows the 6
  // best by default), roc by AUC, hist-clf by each model's own payload F1.
  const cum_styled = $derived(shown(cum_pr.models))
  const box_styled = $derived(order_models(shown(box_data.models), model_mae))
  const roc_styled = $derived(order_models(shown(roc.models), (mdl) => -mdl.auc))
  const hist_clf_styled = $derived(order_models(shown(hist_clf.models), (mdl) => -mdl.f1))
  const rolling_styled = $derived(
    order_models(shown(rolling_mae.models), model_mae).map((mdl, idx) => ({
      ...mdl,
      visible: idx < 6,
    })),
  )

  // each model's cumulative precision/recall line plus a dot marking the end of its
  // stable-prediction ranking (= its total count of predicted-stable materials)
  const cumulative_series = (key: `precision` | `recall`): DataSeries[] =>
    cum_styled.flatMap(({ label, color, x, precision, recall, end }) => [
      {
        x,
        y: key === `precision` ? precision : recall,
        label,
        markers: `line` as const,
        line_style: { stroke: color },
      },
      {
        x: [end[0]],
        y: [end[key === `precision` ? 1 : 2]],
        markers: `points` as const,
        point_style: { fill: color, radius: 2.5 },
      },
    ])
  const stable_count_ref = labeled_vline(cum_pr.n_stable, `Stable Materials`)

  // one shared legend rendered once below both precision+recall plots (a per-plot legend
  // would sit under a single column and squish it). one entry per model line; the
  // unlabeled end-point dot series are intentionally excluded.
  const cum_pr_legend: LegendItem[] = $derived(
    cum_styled.map(({ label, color }, idx) => ({
      label,
      visible: true,
      series_idx: idx * 2,
      display_style: { line_color: color },
    })),
  )

  const roc_series: DataSeries[] = $derived(
    roc_styled.map(({ label, auc, fpr, tpr }) => ({
      x: fpr,
      y: tpr,
      label: `${label} · AUC=${auc.toFixed(2)}`,
      markers: `line` as const,
    })),
  )

  // densify the 3-vertex triangle outline (y = |x|) so spline interpolation of the
  // line and fill region renders straight edges instead of a parabola
  const triangle_x = Array.from({ length: 201 }, (_, idx) => (idx - 100) / 100)
  const triangle = { x: triangle_x, y: triangle_x.map(Math.abs) }

  const rolling_mae_series: DataSeries[] = $derived([
    // 'triangle of peril' outline; the fill region below shades its inside
    {
      ...triangle,
      id: `triangle`,
      label: `MAE > |E<sub>hull dist</sub>|`,
      markers: `line` as const,
      line_style: { stroke: `red`, stroke_width: 1.5 },
    },
    ...rolling_styled.map(({ label, color, y, visible }) => ({
      x: rolling_mae.x,
      y,
      label,
      markers: `line` as const,
      line_style: { stroke: color },
      visible: visible ?? true,
    })),
  ])

  // rolling count of test-set structures per hull-dist bin, drawn as a filled area in a
  // marginal panel above the main plot (sharing its x-axis) rather than overlaid on it
  const density_series: DataSeries[] = [
    {
      ...rolling_mae.density,
      id: `density`,
      label: `Density`,
      markers: `line` as const,
      line_style: { stroke: `rgba(0, 150, 200, 1)` },
    },
  ]
  const density_fill: FillRegion[] = [{
    lower: { type: `constant`, value: 0 },
    upper: { type: `series`, series_id: `density` },
    fill: `rgba(0, 150, 200, 0.15)`,
    curve: `linear`,
  }]
  // shared l/r padding keeps the marginal's x-axis aligned with the main plot's
  const mae_pad = { l: 66, r: 18 }

  // shade the 'triangle of peril': everything above y = |x| (clipped to the y range)
  const peril_region: FillRegion[] = [{
    lower: { type: `series`, series_id: `triangle` },
    upper: { type: `constant`, value: 1 },
    // alpha in the color itself; fill_opacity doesn't reliably apply to plain colors
    fill: `rgba(255, 0, 0, 0.1)`,
    curve: `linear`,
  }]

  const clf_classes = [
    { key: `tp`, label: `True Positive`, color: `lightseagreen` },
    { key: `fn`, label: `False Negative`, color: `orange` },
    { key: `fp`, label: `False Positive`, color: `lightsalmon` },
    { key: `tn`, label: `True Negative`, color: `dodgerblue` },
  ] as const

  // bars must be one bin wide: matterviz defaults bar_width to 0.5 data units, which for
  // ~0.007-wide hull-dist bins overlaps every bar into a solid block (and z-fights the
  // stacked colors). pin it to the bin spacing so adjacent bars touch without overlap.
  const hist_bin_width = hist_clf.bin_centers[1] - hist_clf.bin_centers[0]

  // one shared legend above all 9 histograms (the classes are identical across panels);
  // square swatches match how BarPlot renders bar-series legend entries
  const clf_legend: LegendItem[] = clf_classes.map(({ label, color }, idx) => ({
    label,
    visible: true,
    series_idx: idx,
    display_style: { symbol_type: `Square` as const, symbol_color: color },
  }))
</script>

## Box Plot of Hull Distance Errors

<BoxPlot
series={box_styled.map(({ label, color, quantiles }) => ({ label, color, y: quantiles }))}
whisker_mode="minmax"
y_axis={{ label: `Error in E<sub>hull dist</sub> (eV/atom)`, format: `.3` }}
x_axis={{ tick: { label: { rotation: 90 } } }}
padding={{ b: 144 }}
show_controls={false}
style="height: 600px"
/>

> @label:fig:box-hull-dist-errors Box plot of interquartile ranges (IQR) of hull distance errors for each model. The whiskers extend to the minimum and maximum values. The horizontal line inside the box shows the median. BOWSR has the highest median error, while Voronoi RF has the highest IQR. Note that MEGNet and CGCNN are the only models with a positive median. Their hull distance errors are biased towards more frequently predicting thermodynamic instability, explaining why they are closest to getting the overall number of stable structures in the test set right (see cumulative precision/recall in @fig:cumulative-precision-recall).

## Cumulative Precision and Recall

<div class="fig-grid two-col">
  <figure>
    <figcaption>Precision</figcaption>
    <ScatterPlot
      series={cumulative_series(`precision`)}
      y_axis={{ range: [0, 1] }}
      ref_lines={[stable_count_ref]}
      legend={null}
      style="height: 360px"
    />
  </figure>
  <figure>
    <figcaption>Recall</figcaption>
    <ScatterPlot
      series={cumulative_series(`recall`)}
      y_axis={{ range: [0, 1] }}
      ref_lines={[
        { type: `segment`, p1: [0, 0], p2: [cum_pr.n_stable, 1], style: dashed },
        stable_count_ref,
      ]}
      legend={null}
      style="height: 360px"
    />
  </figure>
</div>

<PlotLegend series_data={cum_pr_legend} draggable={false} filterable={false} style={wide_legend.style} />

> @label:fig:cumulative-precision-recall Model precision and recall for thermodynamic stability classification as a function of number of materials ranked from most to least stable by each model.
> CHGNet initially achieves the highest cumulative precision and recall.
> Simulates materials discovery efforts of different sizes since a typical campaign will rank hypothetical materials by model-predicted hull distance from most to least stable and validate the most stable predictions first.
> A higher fraction of correct stable predictions corresponds to higher precision and fewer stable materials overlooked correspond to higher recall.
> This figure highlights how different models perform better or worse depending on the length of the discovery campaign.
> The UIPs (CHGNet, M3GNet, MACE) are seen to offer significantly improved precision on shorter campaigns of ~20k or less materials validated as they are less prone to false positive predictions among highly stable materials.

## Receiver Operating Characteristic (ROC) Curve

<ScatterPlot
series={roc_series}
legend={wide_legend}
x_axis={{ label: `False Positive Rate`, range: [0, 1.05] }}
y_axis={{ label: `True Positive Rate`, range: [0, 1.05] }}
ref_lines={[{ type: `diagonal`, slope: 1, intercept: 0, style: { ...dashed, color: `gray` } }]}
point_tween={{ duration: 0 }}
line_tween={{ duration: 0 }}
style="height: 700px"
/>

> @label:fig:roc-models Receiver operating characteristic (ROC) curve for each model. TPR/FPR = true/false positive rate. FPR on the $x$ axis is the fraction of unstable structures classified as stable. TPR on the $y$ axis is the fraction of stable structures classified as stable. Dashed diagonal shows a random classifier.

## Rolling MAE vs Hull Distance

<div class="rolling-mae">
  <ScatterPlot
    series={density_series}
    fill_regions={density_fill}
    legend={null}
    x_axis={{ range: [-0.2, 0.2], ticks: [] }}
    y_axis={{ range: [0, 8000], format: `~s` }}
    padding={{ ...mae_pad, t: 5, b: 0 }}
    point_tween={{ duration: 0 }}
    line_tween={{ duration: 0 }}
    style="height: 90px; min-height: 90px"
  />
  <ScatterPlot
    series={rolling_mae_series}
    fill_regions={peril_region}
    legend={wide_legend}
    x_axis={{ label: `E<sub>above MP hull</sub> (eV/atom)`, range: [-0.2, 0.2] }}
    y_axis={{ label: `Rolling MAE (eV/atom)`, range: [0, 0.1] }}
    padding={mae_pad}
    point_tween={{ duration: 0 }}
    line_tween={{ duration: 0 }}
    style="height: 700px"
  />
</div>

> @label:fig:rolling-mae-vs-hull-dist-models Universal potentials are more reliable classifiers because they exit the red triangle earliest.
> These lines show the rolling MAE on the WBM test set as the energy to the convex hull of the MP training set is varied.
> Lower is better.
> Inside the large red 'triangle of peril', models are most likely to misclassify structures.
> As long as a model's rolling MAE remains inside the triangle, its mean error is larger than the distance to the convex hull.
> If the model's error for a given prediction happens to point towards the stability threshold at $E$<sub>above MP hull</sub> = 0, its average error will change the stability classification from true positive/negative to false negative/positive.
> The width of the 'rolling window' box indicates the width over which prediction errors were averaged.

## Distribution of Model-Predicted Hull Distance

<PlotLegend
  series_data={clf_legend}
  layout="horizontal"
  layout_tracks={4}
  draggable={false}
  filterable={false}
  style="margin: 0 auto 1em"
/>

<div class="fig-grid three-col">
  {#each hist_clf_styled as model (model.label)}
  <figure>
    <figcaption>{model.label} · F1={model.f1}</figcaption>
    <BarPlot
      series={clf_classes.map(({ key, label, color }) => ({
        x: hist_clf.bin_centers,
        y: model[key],
        bar_width: hist_bin_width,
        label,
        color,
      }))}
      mode="stacked"
      x_axis={{ range: [-0.4, 0.4] }}
      y_axis={{ range: [0, 18_000], format: `~s` }}
      show_legend={false}
      show_controls={false}
      style="height: 240px"
    />
  </figure>
  {/each}
</div>

> @label:fig:hist-clf-pred-hull-dist-models Distribution of model-predicted hull distance colored by stability classification. Models are sorted from top to bottom by F1 score. The thickness of the red and yellow bands shows how often models misclassify as a function of how far away from the convex hull they place a material. While CHGNet and M3GNet perform almost equally well overall, these plots reveal that they do so via different trade-offs. M3GNet commits fewer false negatives but more false positives predictions compared to CHGNet. In a real discovery campaign, false positives have a higher opportunity cost than false negatives, since they result in wasted DFT relaxations or even synthesis time in the lab. A false negative by contrast is just one missed opportunity out of many. For this reason, models with high true positive rate (TPR) even at the expense of lower true negative rate (TNR) are generally preferred.

<style>
  .fig-grid {
    display: grid;
    gap: 1em 1.5em;
  }
  .fig-grid.two-col {
    grid-template-columns: 1fr 1fr;
  }
  .fig-grid.three-col {
    grid-template-columns: repeat(3, 1fr);
    gap: 0.4em 0.5em;
  }
  .fig-grid figure {
    margin: 0;
  }
  .fig-grid figcaption {
    text-align: center;
    font-size: 0.9em;
    margin-bottom: 0.2em;
  }
  .fig-grid.two-col figcaption {
    margin-top: -0.5em;
  }
  .rolling-mae {
    display: flex;
    flex-direction: column;
    gap: 1.5em;
  }
</style>
