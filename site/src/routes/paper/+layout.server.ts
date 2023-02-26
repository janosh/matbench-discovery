import fs from 'fs'
import yml from 'js-yaml'
import type { LayoutServerLoad } from './$types'

export const load: LayoutServerLoad = async ({ route }) => {
  const data = fs.readFileSync(`src/routes/${route.id}/+page.md`, `utf8`)
  const cff = fs.readFileSync(`../citation.cff`, `utf8`)

  // Count the number of words using a regular expression
  const word_count = data.match(/\b\w+\b/g)?.length ?? null

  return { word_count, ...yml.load(cff) }
}
