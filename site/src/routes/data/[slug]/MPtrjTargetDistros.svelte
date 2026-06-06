<script lang="ts">
  import mp_trj_hists from '$figs/mp-trj-hists.json.gz'
  import { plotly_blue } from '$lib/fig-helpers'
  import { BarPlot } from 'matterviz/plot'
  import MPtrjElemCountsPtable from './MPtrjElemCountsPtable.svelte'
  import MpTrjNSitesHist from './MpTrjNSitesHist.svelte'
  import MPtrjTargetCounts from './mptrj-target-counts.md'

  // per-target presentation: title, x label, and whether the count axis needs
  // arcsinh compression (log-like but keeps a valid 0 baseline for bars)
  const targets = [
    { key: `e-form`, title: `Formation Energy`, x_label: `E<sub>form</sub> (eV/atom)`, arcsinh: false },
    { key: `forces`, title: `Forces`, x_label: `|Forces| (eV/Å)`, arcsinh: true },
    { key: `stresses`, title: `Stresses`, x_label: `1/3 Tr(σ) (eV/Å³)`, arcsinh: true },
    { key: `magmoms`, title: `Magnetic Moments`, x_label: `Magmoms (μ<sub>B</sub>)`, arcsinh: true },
  ] as const
</script>

<MPtrjTargetCounts />

<ul>
  {#each targets as { key, title, x_label, arcsinh } (key)}
    <div>
      <h3>{title}</h3>
      <BarPlot
        series={[{ ...mp_trj_hists[key], color: plotly_blue }]}
        x_axis={{ label: x_label }}
        y_axis={{
          label: `Number of Structures`,
          ...(arcsinh ? { scale_type: `arcsinh` as const } : {}),
        }}
        show_controls={false}
        style="height: 300px; width: 100%; max-width: 700px"
      />
    </div>
  {/each}
  <div>
    <h3>Number of Sites</h3>
    <MpTrjNSitesHist style="height: 300px; width: 100%; max-width: 700px" />
  </div>
</ul>

<h2>Elemental Prevalence</h2>

<MPtrjElemCountsPtable />

<style>
  ul {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
    gap: 1em;
  }
  ul h3 {
    text-transform: capitalize;
    text-align: center;
    margin: 1em auto 0;
  }
</style>
