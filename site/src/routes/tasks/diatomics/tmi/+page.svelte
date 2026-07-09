<script lang="ts">
  import { type InternalPoint, ScatterPlot } from 'matterviz'
  import { element_data, type ChemicalElement } from 'matterviz/element'
  import { ELEM_SYMBOLS } from 'matterviz/labels'
  import { SvelteSet } from 'svelte/reactivity'
  import { make_plot_observer } from '../observe-plot'

  let { data } = $props()
  let magmom_curves = $derived(data?.magmom_curves ?? {})

  // same fixed reference colors as the main diatomics page
  const xc_colors: Record<string, string> = { PBE: `#4c78e0`, r2SCAN: `#f032e6` }
  const elements = element_data as ChemicalElement[]
  const element_by_symbol = new Map(elements.map((element) => [element.symbol, element]))

  type PointMeta = {
    functional: string
    atom_idx: number
    spin_candidate: string | null
  }

  let formulas = $derived(
    ELEM_SYMBOLS.map((symbol) => `${symbol}-${symbol}`).filter(
      (formula) => formula in magmom_curves,
    ),
  )

  const visible_plots = new SvelteSet<string>()
  const observe_plot = make_plot_observer(visible_plots)
  let tooltip_point: InternalPoint | null = $state(null)

  function series_for(formula: string) {
    return Object.entries(magmom_curves[formula] ?? {}).flatMap(([functional, curve]) =>
      [0, 1].map((atom_idx) => {
        const points = curve.distances.flatMap((distance, pt_idx) => {
          const moments = curve.magmoms[pt_idx]
          if (!moments) return []
          return [
            {
              distance,
              magmom: moments[atom_idx],
              spin_candidate: curve.spin_candidates[pt_idx],
            },
          ]
        })
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
    <span><span class="swatch" style="background: {color}"></span>{functional}</span>
  {/each}
</div>

<div class="grid">
  {#each formulas as formula (formula)}
    {@const element_symbol = formula.split(`-`, 1)[0] as ChemicalElement[`symbol`]}
    {@const element = element_by_symbol.get(element_symbol)}
    <div class="plot-shell" {@attach observe_plot(formula)}>
      {#if visible_plots.has(formula)}
        <div class="plot">
          <h3 title={element ? `${element.name} (Z=${element.number})` : formula}>
            {element ? `${element.number} ${formula}` : formula}
          </h3>
          <ScatterPlot
            series={series_for(formula)}
            x_axis={{
              label: `Distance (Å)`,
              format: `.1f`,
              range: [0.2, 6],
              label_shift: { y: -30 },
            }}
            y_axis={{ label: `Magmom (μB)`, format: `.1f` }}
            bind:tooltip_point
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
        </div>
      {:else}
        <div class="plot-placeholder"><h3>{formula}</h3></div>
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
  .legend > span {
    display: inline-flex;
    align-items: center;
    gap: 0.4em;
  }
  .swatch {
    width: 1em;
    height: 0.35em;
    border-radius: 2px;
  }
  .grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(min(100%, 26rem), 1fr));
    gap: 55pt 15pt;
  }
  .plot-shell,
  .plot-placeholder {
    min-height: 300px;
  }
  .plot {
    --scatter-font-size: 14px;
    display: flex;
    flex-direction: column;
    box-sizing: border-box;
    height: 300px;
  }
  .plot-placeholder {
    display: grid;
    place-items: start center;
    color: var(--text-muted, #777);
  }
  h3 {
    align-self: center;
    margin: 0 0 -5pt;
    text-align: center;
    font-size: 0.9em;
  }
</style>
