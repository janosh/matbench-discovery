import { by_date_added_desc, type DiatomicsCurves, MODELS } from '$lib'
import { fetch_diatomics_data } from '$lib/server/diatomics'
import type { PageServerLoad } from './$types'

type PageDiatomicsCurves = {
  distances: number[]
  'homo-nuclear': Record<string, { energies: number[] }>
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
  // Filter models that have diatomics metrics, newest first
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

      const source = diatomics as Record<string, unknown>
      const pred_file =
        typeof source.pred_file === `string` ? source.pred_file : undefined
      const pred_file_url =
        typeof source.pred_file_url === `string` ? source.pred_file_url : undefined

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

  return { diatomic_models, diatomic_curves, errors }
}
