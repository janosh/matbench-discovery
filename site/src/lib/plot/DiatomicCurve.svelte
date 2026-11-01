<script lang="ts" module>
  import type { ChemicalElement } from 'matterviz/element'

  // Plot-title tooltip: full element name, atomic number and periodic-table category
  export const element_tooltip = ({ name, number, category }: ChemicalElement): string =>
    `${name} (Z=${number}) · ${category}`
</script>

<script lang="ts">
  import { ScatterPlot } from 'matterviz/plot'
  import { DEFAULT_CATEGORY_COLORS } from 'matterviz/colors'
  import { element_data } from 'matterviz/element'
  import { format_num } from 'matterviz/labels'
  import { tooltip as add_tooltip } from 'svelte-widgets/attachments'
  import type { HTMLAttributes } from 'svelte/elements'

  let {
    formula,
    curves,
    ...rest
  }: HTMLAttributes<HTMLDivElement> & {
    formula: string
    // magmoms, spin_candidates and bond_forces are only set for DFT references
    curves: {
      model_key: string
      label: string
      distances: number[]
      energies: number[]
      magmoms?: ([number, number] | null)[]
      spin_candidates?: (string | null)[]
      bond_forces?: number[]
      color: string
      line_width?: number
    }[]
  } = $props()

  type PointMeta = {
    magmoms: [number, number] | null
    spin_candidate: string | null
    bond_force: number | null
  }

  // Tooltip number without trailing zeros; anything rounding to zero prints a bare 0
  // (not +0.00 or −0). sign=`+` marks positive values, e.g. FM vs AFM site moments.
  const trim_num = (value: number, decimals: number, sign: `` | `+` = ``): string => {
    const rounded = Number(value.toFixed(decimals)) || 0
    return rounded === 0 ? `0` : format_num(rounded, `${sign}.${decimals}~f`)
  }

  const x_range: [number, number] = [0.2, 6]
  const y_range: [number, number] = [-8, 20]
  const no_tween = { duration: 0 } as const
  const elements = element_data as ChemicalElement[]
  let element_symbol = $derived(formula.split(`-`, 1)[0])
  let element = $derived(
    elements.find((element_entry) => element_entry.symbol === element_symbol),
  )
  let header_text = $derived(element ? `${element.number} ${formula}` : formula)

  let series = $derived(
    curves.map((curve) => {
      // Keep only points within the x range, pairing each distance with its energy
      const x: number[] = []
      const energies: number[] = []
      const metadata: PointMeta[] = []
      curve.distances.forEach((distance, idx) => {
        if (!(distance >= x_range[0] && distance <= x_range[1])) return
        x.push(distance)
        energies.push(curve.energies[idx])
        metadata.push({
          magmoms: curve.magmoms?.[idx] ?? null,
          spin_candidate: curve.spin_candidates?.[idx] ?? null,
          bond_force: curve.bond_forces?.[idx] ?? null,
        })
      })
      const ref_energy = energies.at(-1) ?? 0

      return {
        id: curve.model_key,
        label: curve.label,
        x,
        y: energies.map((energy) => energy - ref_energy),
        metadata,
        markers: `line+points` as const,
        point_style: { fill: curve.color, radius: 2.5, stroke_width: 0 },
        line_style: {
          stroke: curve.color,
          // DFT references (line_width set) get a thicker line to stand out
          ...(curve.line_width ? { stroke_width: curve.line_width } : {}),
        },
      }
    }),
  )
</script>

<div {...rest} class={[`diatomic-plot`, rest.class]}>
  <h3
    class="diatomic-plot-title"
    style:--category-color={element && DEFAULT_CATEGORY_COLORS[element.category]}
    {@attach add_tooltip({
      content: element ? element_tooltip(element) : formula,
      placement: `top`,
    })}
  >
    {header_text}
  </h3>
  <ScatterPlot
    {series}
    x_axis={{
      label: `Distance (Å)`,
      format: `.1f`,
      range: x_range,
      label_shift: { y: 0 },
    }}
    y_axis={{
      label: `Energy (eV)`,
      format: `.2f`,
      range: y_range,
      label_shift: { x: 12 },
    }}
    legend={null}
    point_tween={no_tween}
    line_tween={no_tween}
  >
    {#snippet tooltip({ label, x, y, metadata })}
      {@const { magmoms, spin_candidate, bond_force } = metadata ?? {}}
      <strong>{label ?? ``}</strong><br />
      Distance = {format_num(x, `.3~f`)} Å<br />
      E − E(6 Å) = {trim_num(y, 2)} eV
      {#if bond_force != null}
        <br />Bond force = {trim_num(bond_force, 3)} eV/Å{trim_num(bond_force, 3) === `0`
          ? ``
          : bond_force > 0
            ? ` (repulsive)`
            : ` (attractive)`}
      {/if}
      {#if spin_candidate != null}
        <br />Spin state = {spin_candidate === `afm`
          ? `AFM broken symmetry`
          : `${spin_candidate} μB total (NUPDOWN)`}
      {/if}
      {#if magmoms}
        <br />Site moments = {magmoms
          .map((moment) => trim_num(moment, 2, `+`))
          .join(`, `)} μB
      {/if}
    {/snippet}
  </ScatterPlot>
</div>
