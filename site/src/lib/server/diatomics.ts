// Diatomic curves as served to the diatomics page: one prerendered JSON file per model
// (and per DFT reference) that the page fetches only for selected curves, instead of
// serializing all ~10 MB of curves into the page HTML.
import { by_benchmark_added_desc } from '#lib'
import DFT_REFERENCES from '#lib/diatomics-dft.json.gz'
import { ACTIVE_MODELS, has_diatomics_curves } from '#lib/models.svelte.js'
import { fetch_diatomics_data } from '#lib/server/predictions.js'
import type { DiatomicsCurves } from '#lib/types.js'
import type { PageDiatomicsCurves } from '#routes/benchmarks/diatomics/diatomics.svelte.js'

// newest first: the page's color palette and curve picker follow this order
export const diatomic_models = ACTIVE_MODELS.filter(has_diatomics_curves).toSorted(
  by_benchmark_added_desc,
)
export const reference_names = Object.keys(DFT_REFERENCES)
// every key with a curves/<key>.json file: DFT references plus model keys
export const curve_keys: readonly string[] = [
  ...reference_names,
  ...diatomic_models.map(({ model_key }) => model_key),
]

// plots only need ~6 significant digits (sub-meV near the minimum); full float64
// precision more than doubles the payload
const round = (value: number): number => Number(value.toPrecision(6))

// homonuclear energies only: forces and heteronuclear curves are large and unused
const to_page_curves = (curves: DiatomicsCurves): PageDiatomicsCurves => ({
  distances: curves.distances,
  'homo-nuclear': Object.fromEntries(
    Object.entries(curves[`homo-nuclear`]).map(([formula, curve]) => [
      formula,
      { energies: curve.energies.map(round) },
    ]),
  ),
})

// DFT references carry per-formula distance grids plus the per-point electronic state
// shown in plot tooltips; forces shrink to the bond-axis component on atom 2 (dimers
// lie along x)
const reference_curves = (ref_name: string): PageDiatomicsCurves => ({
  distances: [],
  'homo-nuclear': Object.fromEntries(
    Object.entries(DFT_REFERENCES[ref_name]).map(
      ([formula, { distances, energies, forces, magmoms, spin_candidates }]) => [
        formula,
        {
          distances,
          energies,
          magmoms,
          spin_candidates,
          bond_forces: forces?.map((point_forces) => point_forces[1][0]),
        },
      ],
    ),
  ),
})

// Curves for one model key or DFT reference name. A failed download fails the build
// (like MD) rather than prerendering a page with missing curves.
export async function diatomic_page_curves(key: string): Promise<PageDiatomicsCurves> {
  if (reference_names.includes(key)) return reference_curves(key)
  const model = diatomic_models.find(({ model_key }) => model_key === key)
  if (!model) throw new Error(`No diatomics curves for ${key}`)
  const curves = await fetch_diatomics_data({
    pred_file: model.metrics?.diatomics?.pred_file,
  }).catch((error: unknown) => {
    throw new Error(`Failed to fetch diatomics curves for ${key}`, { cause: error })
  })
  return to_page_curves(curves)
}
