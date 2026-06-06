<script lang="ts">
  import DataReadme from '$data/wbm/readme.md'
  import hist_e_form from '$figs/hist-wbm-e-form-per-atom.json.gz'
  import hist_hull_dist from '$figs/hist-wbm-hull-dist.json.gz'
  import arity_hist from '$figs/mp-vs-mp-trj-vs-wbm-arity-hist.json.gz'
  import sunbursts from '$figs/spacegroup-sunbursts.json.gz'
  import { PtableHeatmap } from '$lib'
  import { dashed, floating_label, labeled_vline } from '$lib/fig-helpers'
  import type { ElementSymbol } from 'matterviz'
  import { ColorScaleSelect, Icon } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'
  import { BarPlot, sunburst_from_labels_parents, Sunburst } from 'matterviz/plot'
  import Select from 'svelte-multiselect'
  import { tooltip } from 'svelte-multiselect/attachments'
  import MPtrjElemCountsPtable from './[slug]/MPtrjElemCountsPtable.svelte'
  import MpTrjNSitesHist from './[slug]/MpTrjNSitesHist.svelte'
  import DataFilesDirectDownload from './data-files-direct-download.md'
  import MpElementalReferenceEnergies from './mp-elemental-reference-energies.md'

  // build the SunburstNode tree from the flat plotly arrays at render time (smaller
  // payload than shipping a pre-nested tree; matterviz handles duplicate ids)
  const build_sunburst = (sb: typeof sunbursts.mp) =>
    sunburst_from_labels_parents(sb.labels, sb.parents, sb.values, { ids: sb.ids })

  const { mean, std } = hist_hull_dist
  const hull_dist_refs = [
    labeled_vline(mean - std, `mean - std = ${(mean - std).toFixed(2)}`),
    labeled_vline(mean, `mean = ${mean.toFixed(2)}`),
    labeled_vline(mean + std, `mean + std = ${(mean + std).toFixed(2)}`),
    floating_label(mean - std, `stable`, `#636efa`),
    floating_label(mean + std, `unstable`, `#ef553b`),
  ]
  const elem_counts = import.meta.glob<Record<ElementSymbol, number>>(
    `./*-element-counts-by-{occurrence,composition}*.json`,
    { eager: true, import: 'default' },
  )

  let log_scale = $state(false) // Log color scale
  let color_scale = $state<D3InterpolateName>(`interpolateViridis`)
  const count_modes = [`occurrence`, `composition`]
  let count_mode = $state(count_modes[0])

  let mp_elem_counts = $derived(
    elem_counts[`./mp-element-counts-by-${count_mode}.json`],
  )
  let mp_trj_elem_counts = $derived(
    elem_counts[`./mp-trj-element-counts-by-${count_mode}.json`],
  )
  let wbm_elem_counts = $derived(
    elem_counts[`./wbm-element-counts-by-${count_mode}.json`],
  )
  $effect.pre(() => {
    if (!mp_elem_counts) throw new Error(`No MP data for count mode ${count_mode}!`)
    if (!mp_trj_elem_counts) throw new Error(`No MPtrj data for count mode ${count_mode}!`)
    if (!wbm_elem_counts) throw new Error(`No WBM data for count mode ${count_mode}!`)
  })

  const capture_state = () => ({ color_scale, log_scale, count_mode })
  export const snapshot = {
    capture: capture_state,
    restore: (
      values: ReturnType<typeof capture_state>,
    ) => ({ color_scale, log_scale, count_mode } = values),
  }
</script>

<DataFilesDirectDownload />

<DataReadme>
  {#snippet hist_e_form_per_atom()}
    <BarPlot
      series={[{ ...hist_e_form, color: `#636efa` }]}
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
    >Count Mode</label>
    <Select
      id="count-mode"
      bind:value={count_mode}
      options={count_modes}
      minSelect={1}
      maxSelect={1}
    >
      {#snippet children({ option })}
        {option}&nbsp;<span
          title="The difference between count modes is best explained by example.
          <code>occurrence</code> mode maps Fe<sub>2</sub>O<sub>3</sub> to <code>{`{Fe: 1, O: 1}`}</code>,
          <code>composition</code> mode maps it to <code>{`{Fe: 2, O: 3}`}</code>."
          {@attach tooltip({ allow_html: true })}
        >
          <Icon icon="Info" style="color: var(--link-color)" />
        </span>
      {/snippet}
    </Select>
    <ColorScaleSelect bind:value={color_scale} selected={[color_scale]} />
    <PtableHeatmap
      heatmap_values={wbm_elem_counts}
      {color_scale}
      colorbar={{
        title: `WBM element counts by ${count_mode}`,
        title_style: `font-size: 1.3em;`,
      }}
      log={log_scale}
    />
  {/snippet}

  {#snippet mp_elements_heatmap()}
    <PtableHeatmap
      heatmap_values={mp_elem_counts}
      {color_scale}
      colorbar={{
        title: `MP element counts by ${count_mode}`,
        title_style: `font-size: 1.3em;`,
      }}
      log={log_scale}
    />
  {/snippet}

  {#snippet mp_trj_elements_heatmap()}
    <MPtrjElemCountsPtable {count_mode} log={log_scale} {color_scale} />
  {/snippet}

  {#snippet hist_wbm_hull_dist()}
    <BarPlot
      series={[
        { ...hist_hull_dist.stable, bar_width: hist_hull_dist.bar_width, color: `#636efa` },
        { ...hist_hull_dist.unstable, bar_width: hist_hull_dist.bar_width, color: `#ef553b` },
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
      <Sunburst data={build_sunburst(sunbursts.mp)} value_mode="total" show_controls={false} style="height: 420px" />
      <Sunburst data={build_sunburst(sunbursts.wbm)} value_mode="total" show_controls={false} style="height: 420px" />
    </div>
  {/snippet}
</DataReadme>

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

<MpTrjNSitesHist style="height: 320px; margin: auto; max-width: 80cqw; padding-right: 2em" />
<p>
  Histogram of number of atoms per structure. The inset shows the same distribution
  log-scaled to visualize the tail of large structures. The green cumulative line in the
  inset shows that 82% have less than 50 sites and 97% of structures in MPtrj have less
  than 100 atoms.
</p>

<style>
  label {
    display: flex;
    gap: 1ex;
    place-content: center;
    align-items: start;
    justify-items: center;
  }
</style>
