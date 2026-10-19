import { by_benchmark_added_desc } from '#lib'
import type { DiatomicsCurves } from '#lib/types.js'
import { ACTIVE_MODELS, has_diatomics_curves } from '#lib/models.svelte.js'
import dft_references from '#lib/diatomics-dft.json.gz'
import { fetch_diatomics_data } from '#lib/server/predictions.js'
import type { PageServerLoad } from './$types'

// VASP PBE/r2SCAN homonuclear references, keyed functional -> formula -> curve. Each
// element has its own distance grid, so curves carry per-formula distances.
const DFT_REFERENCES = dft_references as Record<
  string,
  Record<string, { distances: number[]; energies: number[]; forces?: number[][][] }>
>

type PageDiatomicsCurves = {
  distances: number[]
  // distances is optional per formula (DFT references use per-element grids)
  'homo-nuclear': Record<string, { energies: number[]; distances?: number[] }>
}

const to_page_curves = (curves: DiatomicsCurves): PageDiatomicsCurves => ({
  distances: curves.distances,
  'homo-nuclear': Object.fromEntries(
    Object.entries(curves[`homo-nuclear`]).map(([formula, curve]) => [
      formula,
      { energies: curve.energies },
    ]),
  ),
})

export const load: PageServerLoad = async () => {
  const diatomic_models = ACTIVE_MODELS.filter(has_diatomics_curves).toSorted(
    by_benchmark_added_desc,
  )

  // Fetch data for all models at build time; a failed download fails the build (like
  // MD) rather than prerendering a page with missing curves. Return only the
  // homonuclear energies used by the page; forces and heteronuclear curves are large
  // and unused.
  const diatomic_curves: Record<string, PageDiatomicsCurves> = Object.fromEntries(
    await Promise.all(
      diatomic_models.map(async (model) => {
        const curves = await fetch_diatomics_data({
          pred_file: model.metrics?.diatomics?.pred_file,
        }).catch((error: unknown) => {
          throw new Error(`Failed to fetch diatomics curves for ${model.model_key}`, {
            cause: error,
          })
        })
        return [model.model_key, to_page_curves(curves)] as const
      }),
    ),
  )

  // Add DFT reference curves (PBE, r2SCAN) as selectable pseudo-models on the plots.
  // Keep only distances/energies; forces are used for metrics but not plotted.
  for (const [ref_name, formulas] of Object.entries(DFT_REFERENCES)) {
    diatomic_curves[ref_name] = {
      distances: [],
      'homo-nuclear': Object.fromEntries(
        Object.entries(formulas).map(([formula, { distances, energies }]) => [
          formula,
          { distances, energies },
        ]),
      ),
    }
  }
  const reference_names = Object.keys(DFT_REFERENCES)

  // Model keys only: the page resolves them from MODELS so it sorts by live (re-weighted)
  // scores rather than stale serialized copies.
  const diatomic_model_keys = diatomic_models.map(({ model_key }) => model_key)
  return { diatomic_model_keys, diatomic_curves, reference_names }
}
