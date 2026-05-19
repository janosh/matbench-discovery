import { type DiatomicsCurves, MODELS } from '$lib'
import { fetch_diatomics_data } from '$lib/server/diatomics'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = async () => {
  // Filter models that have diatomics metrics
  const diatomic_models = MODELS.filter(
    (model) => model.metrics?.diatomics && typeof model.metrics.diatomics === `object`,
  ).toSorted(
    (m1, m2) => new Date(m2.date_added).getTime() - new Date(m1.date_added).getTime(), // Sort by date added, newest first
  )

  // Fetch data for all models at build time
  const diatomic_curves: Record<string, DiatomicsCurves> = {}
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
        diatomic_curves[model.model_name] = await fetch_diatomics_data({
          pred_file,
          pred_file_url,
        })
      } catch (error) {
        console.error(`Failed to fetch data for ${model.model_name}:`, error)
        errors[model.model_name] = error instanceof Error ? error.message : String(error)
      }
    }),
  )

  return { diatomic_models, diatomic_curves, errors }
}
