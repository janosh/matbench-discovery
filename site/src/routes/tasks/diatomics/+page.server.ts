import { MODEL_METADATA, type DiatomicsCurves } from '$lib'
import type { PageServerLoad } from './$types'

async function fetch_diatomics(file_url: string): Promise<DiatomicsCurves> {
  const response = await fetch(file_url)

  if (!response.ok) {
    throw new Error(
      `request for ${file_url} failed: ${response.status} ${response.statusText}`,
    )
  }

  // Get response as a stream and decompress it
  const ds = new DecompressionStream(`gzip`)
  const decompressed_stream = response.body?.pipeThrough(ds)
  if (!decompressed_stream) {
    throw new Error(`Failed to decompress response`)
  }

  // Convert the stream to text
  const decompressed_response = new Response(decompressed_stream)
  const text = await decompressed_response.text()
  return JSON.parse(text)
}

export const load: PageServerLoad = async () => {
  // Filter models that have diatomics metrics
  const diatomic_models = MODEL_METADATA.filter(
    (model) => model.metrics?.diatomics && typeof model.metrics.diatomics === `object`,
  ).sort(
    (m1, m2) => new Date(m2.date_added).getTime() - new Date(m1.date_added).getTime(), // Sort by date added, newest first
  )

  // Fetch data for all models at build time
  const diatomic_curves: Record<string, DiatomicsCurves> = {}
  const errors: Record<string, string> = {}

  await Promise.all(
    diatomic_models.map(async (model) => {
      const { diatomics } = model.metrics
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
      } catch (err) {
        console.error(`Failed to fetch data for ${model.model_name}:`, err)
        errors[model.model_name] = err instanceof Error ? err.message : String(err)
      }
    }),
  )

  return { diatomic_models, diatomic_curves, errors }
}
