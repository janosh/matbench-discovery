<script lang="ts">
  import { browser } from '$app/environment'
  import DataReadme from '$data/wbm/readme.md'
  import WbmFormEnergyHist from '$figs/hist-wbm-e-form-per-atom.svelte'
  import HistWbmHullDist from '$figs/hist-wbm-hull-dist.svelte'
  import MPtrjNSitesHist from '$figs/mp-trj-n-sites-hist.svelte'
  import MPvsMPtrjVsWBMArityHist from '$figs/mp-vs-mp-trj-vs-wbm-arity-hist.svelte'
  import SpacegroupSunburstMp from '$figs/spacegroup-sunburst-mp.svelte'
  import SpacegroupSunburstWbm from '$figs/spacegroup-sunburst-wbm.svelte'
  import { Icon, PtableHeatmap } from '$lib'
  import type { ElementSymbol } from 'matterviz'
  import { ColorScaleSelect } from 'matterviz'
  import type { D3InterpolateName } from 'matterviz/colors'
  import Select from 'svelte-multiselect'
  import { tooltip } from 'svelte-multiselect/attachments'
  import MPtrjElemCountsPtable from './[slug]/MPtrjElemCountsPtable.svelte'
  import DataFilesDirectDownload from './data-files-direct-download.md'
  import MpElementalReferenceEnergies from './mp-elemental-reference-energies.md'

  const elem_counts: Record<string, Record<ElementSymbol, number>> = import.meta.glob(
    `./*-element-counts-by-{occurrence,composition}*.json`,
    { eager: true, import: `default` },
  )

  let log_scale = $state(false) // log color scale
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
    if (!mp_elem_counts) throw `No MP data for count mode ${count_mode}!`
    if (!mp_trj_elem_counts) throw `No MPtrj data for count mode ${count_mode}!`
    if (!wbm_elem_counts) throw `No WBM data for count mode ${count_mode}!`
  })

  export const snapshot = {
    capture: () => ({ color_scale, log_scale, count_mode }),
    restore: (values) => ({ color_scale, log_scale, count_mode } = values),
  }
</script>

<DataFilesDirectDownload />

<DataReadme>
  {#snippet hist_e_form_per_atom()}
    {#if browser}
      <WbmFormEnergyHist />
    {/if}
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
          {@attach tooltip()}
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
    {#if browser}
      <HistWbmHullDist />
    {/if}
  {/snippet}

  {#snippet spacegroup_sunbursts()}
    {#if browser}
      <div style="display: flex; gap: 2em; place-content: center">
        <SpacegroupSunburstMp />
        <SpacegroupSunburstWbm />
      </div>
    {/if}
  {/snippet}
</DataReadme>

<MpElementalReferenceEnergies />

{#if browser}
  <MPvsMPtrjVsWBMArityHist style="margin: auto; max-width: 60cqw; padding-right: 2em" />
{/if}
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

<MPtrjNSitesHist style="margin: auto; max-width: 80cqw; padding-right: 2em" />
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
