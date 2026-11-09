import { curve_keys, diatomic_page_curves } from '#lib/server/diatomics.js'
import { error } from '@sveltejs/kit'
import type { EntryGenerator, RequestHandler } from './$types'

export const prerender = true

export const entries: EntryGenerator = () => curve_keys.map((key) => ({ key }))

export const GET: RequestHandler = async ({ params }) => {
  if (!curve_keys.includes(params.key)) error(404, `No diatomics curves for ${params.key}`)
  return Response.json(await diatomic_page_curves(params.key))
}
