import { by_date_added_desc, type DiatomicsCurves, MODELS } from '$lib'
import dft_references from '$lib/diatomics-dft.json.gz'
import { fetch_diatomics_data } from '$lib/server/diatomics'
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
  const diatomic_models = MODELS.filter(
    (model) => model.metrics?.diatomics && typeof model.metrics.diatomics === `object`,
  ).toSorted(by_date_added_desc)

  // Fetch data for all models at build time. Return only the homonuclear
  // energies used by the page; forces and heteronuclear curves are large and unused.
  const diatomic_curves: Record<string, PageDiatomicsCurves> = {}
  const errors: Record<string, string> = {}

  await Promise.all(
    diatomic_models.map(async (model) => {
      const diatomics = model.metrics?.diatomics
      if (typeof diatomics !== `object` || diatomics === null) return

      const pred_file =
        typeof diatomics.pred_file === `string` ? diatomics.pred_file : undefined
      const pred_file_url =
        typeof diatomics.pred_file_url === `string` ? diatomics.pred_file_url : undefined

      if (!pred_file && !pred_file_url) {
        errors[model.model_name] = `No prediction file path or URL`
        return
      }

      try {
        const curves = await fetch_diatomics_data({ pred_file, pred_file_url })
        diatomic_curves[model.model_name] = to_page_curves(curves)
      } catch (error) {
        console.error(`Failed to fetch data for ${model.model_name}:`, error)
        errors[model.model_name] = error instanceof Error ? error.message : String(error)
      }
    }),
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

  return { diatomic_models, diatomic_curves, errors, reference_names }
}
