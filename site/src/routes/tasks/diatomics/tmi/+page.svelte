<script lang="ts">
  import SelectToggle from '$lib/SelectToggle.svelte'
  import { bind_url_params, valid_query_param } from '$lib/url-state.svelte'
  import { ScatterPlot } from 'matterviz'
  import type { ChemicalElement } from 'matterviz/element'
  import { ELEM_SYMBOLS } from 'matterviz/labels'
  import { SvelteSet } from 'svelte/reactivity'
  import {
    element_by_symbol,
    element_group_keys,
    element_groups,
  } from '../element-groups'
  import { make_plot_observer } from '../observe-plot'

  type MagmomCurve = {
    distances: number[]
    magmoms: ([number, number] | null)[]
    spin_candidates: (string | null)[]
  }

  let {
    data,
  }: {
    data: { magmom_curves: Record<string, Record<string, MagmomCurve>> }
  } = $props()
  let magmom_curves = $derived(data?.magmom_curves ?? {})

  // same fixed reference colors as the main diatomics page
  const xc_colors: Record<string, string> = { PBE: `#4c78e0`, r2SCAN: `#f032e6` }
  const visible_functionals = new SvelteSet(Object.keys(xc_colors))

  type PointMeta = {
    functional: string
    atom_idx: number
    spin_candidate: string | null
  }

  let selected_element_group = $state(`all`)
  let selected_group = $derived(
    element_groups.find((group) => group.value === selected_element_group) ??
      element_groups[0],
  )
  let formulas = $derived(
    ELEM_SYMBOLS.filter((symbol) => {
      const formula = `${symbol}-${symbol}`
      const element = element_by_symbol.get(symbol)
      return (
        formula in magmom_curves && Boolean(element && selected_group.includes(element))
      )
    }).map((symbol) => `${symbol}-${symbol}`),
  )

  const visible_plots = new SvelteSet<string>()
  const observe_plot = make_plot_observer(visible_plots)

  bind_url_params(
    (params) => {
      selected_element_group = valid_query_param(
        params,
        `elements`,
        `all`,
        element_group_keys,
      )
    },
    () => [[`elements`, selected_element_group, `all`]],
  )

  function series_for(formula: string) {
    return Object.entries(magmom_curves[formula] ?? {})
      .filter(([functional]) => visible_functionals.has(functional))
      .flatMap(([functional, curve]) =>
        [0, 1].map((atom_idx) => {
          const points = curve.distances.map((distance, point_idx) => ({
            distance,
            magmom: curve.magmoms[point_idx]?.[atom_idx] ?? Number.NaN,
            spin_candidate: curve.spin_candidates[point_idx],
          }))
          const color = xc_colors[functional] ?? `gray`
          return {
            id: `${formula}-${functional}-atom${atom_idx + 1}`,
            label: `${functional} atom ${atom_idx + 1}`,
            x: points.map((point) => point.distance),
            y: points.map((point) => point.magmom),
            markers: `line+points` as const,
            metadata: points.map(
              (point): PointMeta => ({
                functional,
                atom_idx: atom_idx + 1,
                spin_candidate: point.spin_candidate,
              }),
            ),
            point_style: {
              fill: color,
              radius: 1.5,
              stroke_width: 0,
              ...(atom_idx === 1 ? { fill_opacity: 0.5 } : {}),
            },
            line_style: {
              stroke: color,
              // atom 2 dashed so overlapping FM branches (m1 == m2) stay legible
              ...(atom_idx === 1 ? { line_dash: `4 3` } : {}),
            },
          }
        }),
      )
  }
</script>

<h1>Diatomics TMI: DFT Reference Spin States</h1>

<p>
  Site-projected magnetic moment of each atom (VASP <code>LORBIT=11</code>, in μB) versus
  dimer separation for the
  <a href="/tasks/diatomics">bundled PBE/r2SCAN DFT reference curves</a>. Solid lines show
  atom 1, dashed lines atom 2 (they overlap for ferromagnetically coupled branches and
  mirror each other for broken-symmetry AFM ones). The per-distance spin-ladder winner is
  shown in the tooltip. Discontinuities in these atom-wise moments flag SCF spin-state
  hops that the total energy curve can hide &mdash; a debugging view suggested by
  <a href="https://cbe.princeton.edu/people/andrew-rosen">Andrew Rosen</a>. Full
  provenance in the
  <a
    href="https://github.com/janosh/matbench-discovery/blob/main/site/src/lib/diatomics-dft.readme.md"
    >DFT reference readme</a
  >.
</p>

<div class="legend">
  {#each Object.entries(xc_colors) as [functional, color] (functional)}
    {@const visible = visible_functionals.has(functional)}
    <button
      class:visible
      aria-pressed={visible}
      onclick={() =>
        visible
          ? visible_functionals.delete(functional)
          : visible_functionals.add(functional)}
    >
      <span class="swatch" style="background: {color}"></span>{functional}
    </button>
  {/each}
</div>

<div class="controls">
  <div aria-label="Element group filter">
    <SelectToggle options={element_groups} bind:selected={selected_element_group} />
  </div>
</div>

<div class="diatomics-grid bleed-1400">
  {#each formulas as formula (formula)}
    {@const element_symbol = formula.split(`-`, 1)[0] as ChemicalElement[`symbol`]}
    {@const element = element_by_symbol.get(element_symbol)}
    <div
      class="diatomic-plot-shell"
      class:diatomic-plot={visible_plots.has(formula)}
      class:diatomic-plot-placeholder={!visible_plots.has(formula)}
      {@attach observe_plot(formula)}
    >
      {#if visible_plots.has(formula)}
        <h2
          class="diatomic-plot-title toc-exclude"
          title={element ? `${element.name} (Z=${element.number})` : formula}
        >
          {element ? `${element.number} ${formula}` : formula}
        </h2>
        <ScatterPlot
          series={series_for(formula)}
          x_axis={{
            label: `Distance (Å)`,
            format: `.1f`,
            range: [0.2, 6],
            label_shift: { y: -30 },
          }}
          y_axis={{ label: `Magmom (μB)`, format: `.1f` }}
          legend={null}
          point_tween={{ duration: 0 }}
          line_tween={{ duration: 0 }}
        >
          {#snippet tooltip({ x_formatted, y_formatted, metadata })}
            <strong>{metadata?.functional} atom {metadata?.atom_idx}</strong><br />
            Distance = {x_formatted} Å<br />
            Magmom = {y_formatted} μB<br />
            {#if metadata?.spin_candidate != null}
              Spin candidate = {metadata.spin_candidate === `afm`
                ? `AFM`
                : `NUPDOWN=${metadata.spin_candidate}`}
            {/if}
          {/snippet}
        </ScatterPlot>
      {:else}
        <h2 class="diatomic-plot-title toc-exclude">{formula}</h2>
      {/if}
    </div>
  {/each}
</div>

<style>
  .legend {
    display: flex;
    gap: 1.5em;
    justify-content: center;
    margin: 1em 0;
  }
  .legend button {
    display: inline-flex;
    align-items: center;
    gap: 0.4em;
    padding: 2pt 5pt;
    background: transparent;
    opacity: 0.45;
  }
  .legend button.visible {
    opacity: 1;
  }
  .legend button:hover {
    background: var(--btn-bg-hover);
  }
  .controls {
    margin-block: 1em;
  }
  .swatch {
    width: 1em;
    height: 0.35em;
    border-radius: 2px;
  }
  .diatomic-plot {
    height: 300px;
  }
</style>
