<script lang="ts">
  import DataReadme from '$data/wbm/readme.md'
  import hist_e_form from '$figs/hist-wbm-e-form-per-atom.json.gz'
  import hist_hull_dist from '$figs/hist-wbm-hull-dist.json.gz'
  import arity_hist from '$figs/mp-vs-mp-trj-vs-wbm-arity-hist.json.gz'
  import sunbursts from '$figs/spacegroup-sunbursts.json.gz'
  import PtableHeatmap from '$lib/PtableHeatmap.svelte'
  import {
    dashed,
    floating_label,
    labeled_vline,
    series_blue,
    series_red,
  } from '$lib/fig-helpers'
  import type { ElementSymbol } from 'matterviz'
  import {
    ColorScaleSelect,
    BarPlot,
    sunburst_from_labels_parents,
    Sunburst,
  } from 'matterviz/plot'
  import { format_num } from 'matterviz/labels'
  import type { D3InterpolateName } from 'matterviz/colors'
  import { Icon, MultiSelect, Popover } from 'svelte-widgets'
  import { Info } from 'svelte-widgets/icons'
  import { bind_url_params, url_color_scale } from '$lib/url-state.svelte'
  import {
    bool_from_param,
    bool_url_entry,
    valid_query_param,
  } from 'svelte-widgets/url-params'
  import MPtrjElemCountsPtable from './MPtrjElemCountsPtable.svelte'
  import MpTrjNSitesHist from './MpTrjNSitesHist.svelte'
  import MpElementalReferenceEnergies from '../mp-elemental-reference-energies.md'

  const { mean, std } = hist_hull_dist
  const hull_dist_refs = [
    labeled_vline(mean - std, `mean - std = ${format_num(mean - std, `.2f`)}`),
    labeled_vline(mean, `mean = ${format_num(mean, `.2f`)}`),
    labeled_vline(mean + std, `mean + std = ${format_num(mean + std, `.2f`)}`),
    floating_label(mean - std, `stable`, series_blue),
    floating_label(mean + std, `unstable`, series_red),
  ]
  const elem_counts = import.meta.glob<Record<ElementSymbol, number>>(
    `../*-element-counts-by-{occurrence,composition}*.json`,
    { eager: true, import: 'default' },
  )

  let log_scale = $state(false) // Log color scale
  let color_scale = $state<D3InterpolateName>(url_color_scale.default)
  const count_modes = [`occurrence`, `composition`]
  let count_mode = $state(count_modes[0])

  const read_url_params = (params: URLSearchParams) => {
    count_mode = valid_query_param(
      params,
      `count_mode`,
      count_modes[0],
      new Set(count_modes),
    )
    log_scale = bool_from_param(params, `log`)
    color_scale = url_color_scale.read(params)
  }
  bind_url_params(read_url_params, () => [
    [`count_mode`, count_mode, count_modes[0]],
    bool_url_entry(`log`, log_scale),
    url_color_scale.entry(color_scale),
  ])

  $effect.pre(() => {
    for (const dataset of [`mp`, `mp-trj`, `wbm`]) {
      if (!elem_counts[`../${dataset}-element-counts-by-${count_mode}.json`])
        throw new Error(`No ${dataset} data for count mode ${count_mode}!`)
    }
  })
</script>

{#snippet elements_heatmap(dataset: `mp` | `wbm`)}
  <PtableHeatmap
    heatmap_values={elem_counts[`../${dataset}-element-counts-by-${count_mode}.json`]}
    {color_scale}
    colorbar={{
      title: `${dataset.toUpperCase()} element counts by ${count_mode}`,
      title_style: `font-size: 1.3em;`,
    }}
    bind:log={log_scale}
  />
{/snippet}

<DataReadme>
  {#snippet hist_e_form_per_atom()}
    <BarPlot
      series={[{ ...hist_e_form, color: series_blue }]}
      x_axis={{ label: `WBM uncorrected formation energy (eV/atom)` }}
      y_axis={{ label: `Number of Structures`, scale_type: `arcsinh` }}
      ref_lines={[
        { type: `vertical`, x: -5, style: dashed },
        { type: `vertical`, x: 5, style: dashed },
      ]}
      show_controls={false}
      style="height: 320px"
    />
  {/snippet}

  {#snippet wbm_elements_heatmap()}
    <label
      for="count-mode"
      style="display: inline-block; transform: translate(10cqw, 5ex)"
    >
      Count Mode
    </label>
    <MultiSelect
      id="count-mode"
      bind:value={count_mode}
      options={count_modes}
      mode="single"
      min_select={1}
    >
      {#snippet children({ option })}
        {option}&nbsp;<Popover
          trigger_mode="hover"
          trap_focus={false}
          aria-label="Count modes"
        >
          {#snippet trigger(trigger_props)}
            <span role="button" tabindex="0" {...trigger_props}
              ><Icon icon={Info} style="color: var(--link-color)" /></span
            >
          {/snippet}
          The difference between count modes is best explained by example.
          <code>occurrence</code> mode maps Fe<sub>2</sub>O<sub>3</sub> to
          <code>{'{Fe: 1, O: 1}'}</code>,
          <code>composition</code> mode maps it to <code>{'{Fe: 2, O: 3}'}</code>.
        </Popover>
      {/snippet}
    </MultiSelect>
    <ColorScaleSelect bind:value={color_scale} aria-label="Color scale" />
    {@render elements_heatmap(`wbm`)}
  {/snippet}

  {#snippet mp_elements_heatmap()}
    {@render elements_heatmap(`mp`)}
  {/snippet}

  {#snippet mp_trj_elements_heatmap()}
    <MPtrjElemCountsPtable {count_mode} bind:log={log_scale} {color_scale} />
  {/snippet}

  {#snippet hist_wbm_hull_dist()}
    <BarPlot
      series={[
        {
          ...hist_hull_dist.stable,
          bar_width: hist_hull_dist.bar_width,
          color: series_blue,
        },
        {
          ...hist_hull_dist.unstable,
          bar_width: hist_hull_dist.bar_width,
          color: series_red,
        },
      ]}
      x_axis={{ label: `WBM energy above MP convex hull (eV/atom)` }}
      y_axis={{ label: `Number of Structures`, format: `~s` }}
      ref_lines={hull_dist_refs}
      show_controls={false}
      style="height: 320px"
    />
  {/snippet}

  {#snippet spacegroup_sunbursts()}
    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 2em">
      {#each [sunbursts.mp, sunbursts.wbm] as { labels, parents, values, ids }}
        <Sunburst
          data={sunburst_from_labels_parents(labels, parents, values, { ids })}
          value_mode="total"
          show_controls={false}
          style="height: 420px"
        />
      {/each}
    </div>
  {/snippet}
</DataReadme>

<p>
  <a href="/data/tmi"
    >Explore WBM element coverage by composition arity and substitution batch.</a
  >
</p>

<MpElementalReferenceEnergies />

<BarPlot
  series={arity_hist.datasets}
  mode="grouped"
  x_axis={{ label: `Number of Elements in Formula` }}
  y_axis={{ label: `Fraction of Structures in Dataset` }}
  show_legend
  show_controls={false}
  style="height: 320px; margin: auto; max-width: 60cqw; padding-right: 2em"
/>
<p>
  Distribution of unique elements per structure in MP, MPtrj and WBM. The bar heights are
  normalized by the total number of structures in each data set. WBM is dominated by
  ternary phases making up 74% of the data set followed by about 13% each of binaries and
  quaternaries. MP has a more even distribution, in particular with more than double the
  relative share of quaternary phases and a significant number of quinternaries which are
  almost absent from WBM. Not shown in this plot for visual clarity are 3% of MP
  structures containing more than 5 elements (up to 9). We also include MPtrj in this plot
  to show a slight drop in relative abundance of quinternary and higher phases vs MP
  ground states.
</p>

<MpTrjNSitesHist
  style="height: 320px; margin: auto; max-width: 80cqw; padding-right: 2em"
/>
<p>
  Histogram of number of atoms per structure. The inset shows the same distribution
  log-scaled to visualize the tail of large structures. The green cumulative line in the
  inset shows that 82% have less than 50 sites and 97% of structures in MPtrj have less
  than 100 atoms.
</p>
