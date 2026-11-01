import { by_benchmark_added_desc } from '#lib'
import type { DiatomicsCurves } from '#lib/types.js'
import { ACTIVE_MODELS, has_diatomics_curves } from '#lib/models.svelte.js'
import DFT_REFERENCES from '#lib/diatomics-dft.json.gz'
import { fetch_diatomics_data } from '#lib/server/predictions.js'
import type { PageServerLoad } from './$types'

// DFT references add per-formula distance grids plus the per-point electronic state
// shown in plot tooltips; model predictions carry energies only
type PageDiatomicCurve = {
  energies: number[]
  distances?: number[]
  magmoms?: ([number, number] | null)[] // site-projected moments of the two atoms (μB)
  spin_candidates?: (string | null)[] // winning NUPDOWN (as string) or `afm`
  bond_forces?: number[] // force on atom 2 along the bond (eV/Å), > 0 = repulsive
}
type PageDiatomicsCurves = {
  distances: number[]
  'homo-nuclear': Record<string, PageDiatomicCurve>
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
  // Each element has its own distance grid, so curves carry per-formula distances.
  // Forces shrink to the bond-axis component on atom 2 (dimers lie along x).
  for (const [ref_name, formulas] of Object.entries(DFT_REFERENCES)) {
    diatomic_curves[ref_name] = {
      distances: [],
      'homo-nuclear': Object.fromEntries(
        Object.entries(formulas).map(
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
    }
  }
  const reference_names = Object.keys(DFT_REFERENCES)

  // Model keys only: the page resolves them from MODELS so it sorts by live (re-weighted)
  // scores rather than stale serialized copies.
  const diatomic_model_keys = diatomic_models.map(({ model_key }) => model_key)
  return { diatomic_model_keys, diatomic_curves, reference_names }
}
