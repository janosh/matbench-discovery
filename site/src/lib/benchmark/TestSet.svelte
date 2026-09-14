<script lang="ts">
  import { benchmarks, type BenchmarkTask } from './data'
  import summary from '$routes/data/benchmark-reference-summary.json'
  import hull_dist from '$figs/hist-wbm-hull-dist.json.gz'
  import spacegroups from '$figs/spacegroup-sunbursts.json.gz'
  import { format_num } from 'matterviz/labels'
  import ElementHeatmap from './ElementHeatmap.svelte'
  import { BarPlot, Histogram } from 'matterviz/plot'
  import { Icon } from 'svelte-widgets'
  import type { Snippet } from 'svelte'

  let { task, children }: { task: BenchmarkTask; children?: Snippet } = $props()
  const dataset = $derived(benchmarks[task].dataset)
  const crystal_systems = spacegroups.wbm.labels.flatMap((label, idx) =>
    spacegroups.wbm.parents[idx] === ``
      ? [{ label, count: spacegroups.wbm.values[idx] }]
      : [],
  )
  const structure_distributions = [
    {
      title: `WBM structure sizes`,
      axis_label: `Atoms per unit cell`,
      scale_type: `log` as const,
      ...summary.wbm_n_sites,
    },
    {
      title: `WBM initial crystal systems`,
      axis_label: `Initial crystal system`,
      scale_type: `linear` as const,
      x: crystal_systems.map(({ label }) => label),
      y: crystal_systems.map(({ count }) => count),
    },
  ]
</script>

<section aria-labelledby="test-set">
  <h2 id="test-set">Test set: {dataset.name}</h2>
  <div class="test-set-summary">
    <div>
      <p class="coverage">{dataset.coverage}</p>
      {#if dataset.credit}
        <p class="credit">Data by <a href={dataset.credit[1]}>{dataset.credit[0]}</a></p>
      {/if}
      <p>{dataset.description}</p>
      <p class="availability">{dataset.availability}</p>
      {#if dataset.id === `wbm`}
        <p>
          Shared by <a href="/benchmarks/discovery#test-set">discovery</a> and
          <a href="/benchmarks/geo-opt#test-set">geometry optimization</a>. These
          summaries describe the full WBM set.
        </p>
      {/if}
      <div class="links">
        {#each dataset.links as [label, href, icon] (label)}
          <a {href}><Icon {icon} /> {label}</a>
        {/each}
      </div>
    </div>
    <figure aria-label="{dataset.name} element occurrences">
      <ElementHeatmap {dataset} />
      <figcaption>
        Each element is counted once per {dataset.count_unit.toLowerCase().slice(0, -1)}.
      </figcaption>
    </figure>
  </div>

  {#if task === `discovery`}
    <figure aria-label="WBM reference hull-distance distribution">
      <BarPlot
        series={[
          {
            ...hull_dist.stable,
            label: `Stable`,
            color: `#4c78a8`,
            bar_width: hull_dist.bar_width,
          },
          {
            ...hull_dist.unstable,
            label: `Unstable`,
            color: `#e45756`,
            bar_width: hull_dist.bar_width,
          },
        ]}
        x_axis={{ label: `DFT energy above MP convex hull (eV/atom)` }}
        y_axis={{ label: `Structures`, format: `~s` }}
        show_controls={false}
        style="height: 300px"
      />
      <figcaption>
        DFT reference energies within two standard deviations of the mean. Negative values
        lie below the fixed MP hull.
      </figcaption>
    </figure>
  {:else if task === `geo-opt`}
    <div class="distributions">
      {#each structure_distributions as { title, axis_label, scale_type, ...series } (title)}
        <figure aria-label={title}>
          <BarPlot
            series={[series]}
            x_axis={{ label: axis_label }}
            y_axis={{ label: `Structures`, format: `~s`, scale_type }}
            show_controls={false}
            style="height: 300px"
          />
        </figure>
      {/each}
    </div>
  {:else if task === `phonons`}
    <figure aria-label="PhononDB reference conductivity distribution">
      <Histogram
        series={[{ values: summary.phonondb_kappa }]}
        bins={30}
        x_axis={{ label: `DFT thermal conductivity at 300 K (W/m/K)`, scale_type: `log` }}
        y_axis={{ label: `Structures` }}
        show_controls={false}
        style="height: 300px"
      />
      <figcaption>
        Direction-averaged reference conductivity, with logarithmic bins to show the range
        across materials.
      </figcaption>
    </figure>
  {:else if task === `md`}
    <!-- svelte-ignore a11y_no_noninteractive_tabindex (Keyboard users must be able to scroll the table.) -->
    <div
      class="systems"
      tabindex="0"
      role="region"
      aria-label="DynaMat reference systems"
    >
      <table>
        <thead
          ><tr
            ><th>System</th><th>Composition</th><th>Atoms</th><th>Temperature (K)</th><th
              >Reference duration (ps)</th
            ></tr
          ></thead
        >
        <tbody>
          {#each summary.md_systems as system (system.name)}
            <tr>
              <td title={system.name}>{system.name.split(`_`)[0].replace(/^bulk/, ``)}</td
              >
              <td>{system.formula}</td><td>{system.n_atoms}</td><td
                >{format_num(system.temperature)}</td
              ><td>{format_num(system.duration_ps, `.2f`)}</td>
            </tr>
          {/each}
        </tbody>
      </table>
    </div>
  {/if}
  {@render children?.()}
</section>

<style>
  section {
    margin-block: 2.5em;
  }
  .test-set-summary,
  .distributions {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    align-items: start;
    gap: 2em;
  }
  .coverage {
    font-weight: 600;
  }
  .credit,
  .availability,
  figcaption {
    color: var(--text-secondary);
    font-size: 0.9em;
  }
  figure {
    min-width: 0;
    margin: 1em 0;
  }
  figcaption {
    margin-top: 0.5em;
    text-align: center;
  }
  .links {
    display: flex;
    flex-wrap: wrap;
    gap: 0.5em 1em;
  }
  .links a {
    display: inline-flex;
    align-items: center;
    gap: 0.3em;
  }
  .systems {
    overflow-x: auto;
  }
  table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.9em;
  }
  th,
  td {
    padding: 0.4em 0.6em;
    text-align: left;
  }
  th {
    white-space: nowrap;
  }
  tbody tr:nth-child(odd) {
    background: var(--card-bg);
  }
  @media (max-width: 700px) {
    .test-set-summary,
    .distributions {
      grid-template-columns: 1fr;
      gap: 0.5em;
    }
  }
</style>
