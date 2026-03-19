// deno-lint-ignore-file no-await-in-loop
import { type DiatomicsCurves, MODELS } from '$lib'
import type { PageServerLoad } from './$types'

async function fetch_diatomics(file_url: string, retries = 3): Promise<DiatomicsCurves> {
  let last_err: unknown
  for (let attempt = 0; attempt < retries; attempt++) {
    try {
      const response = await fetch(file_url)
      if (!response.ok) throw new Error(`${response.status} ${response.statusText}`)

      const decompressed = response.body?.pipeThrough(new DecompressionStream(`gzip`))
      if (!decompressed) throw new Error(`Failed to decompress response`)
      return JSON.parse(await new Response(decompressed).text())
    } catch (error) {
      last_err = error
      if (attempt < retries - 1) {
        // exponential backoff: 1s, 2s, 4s
        await new Promise((resolve) => setTimeout(resolve, 1000 * 2 ** attempt))
      }
    }
  }
  throw new Error(`${file_url} failed after ${retries} attempts: ${String(last_err)}`)
}

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
      if (typeof diatomics !== `object` || !diatomics?.pred_file_url) {
        errors[model.model_name] = `No prediction file URL`
        return
      }

      try {
        const { pred_file_url } = diatomics
        if (typeof pred_file_url !== `string`) {
          errors[model.model_name] = `Invalid prediction file URL`
          return
        }
        diatomic_curves[model.model_name] = await fetch_diatomics(pred_file_url)
      } catch (error) {
        console.error(`Failed to fetch data for ${model.model_name}:`, error)
        errors[model.model_name] = error instanceof Error ? error.message : String(error)
      }
    }),
  )

  return { diatomic_models, diatomic_curves, errors }
}
