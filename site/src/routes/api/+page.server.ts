import { building } from '$app/env'
import { existsSync, readdirSync, readFileSync } from 'node:fs'
import { default_highlighter } from 'svelte-widgets/highlight'
import { assert_ok, create_markdown } from 'svelte-widgets/markdown'

// Generated from docstrings by scripts/make_api_docs.py (optional in dev, required to build).
// Read from disk and rendered as plain Markdown instead of compiled as Svelte, so
// docstring braces like `({} if absent)` stay literal.
const modules_dir = `src/routes/api/modules`

export const load = async () => {
  const markdown = create_markdown({ highlight: default_highlighter.highlight })
  const files = existsSync(modules_dir) ? readdirSync(modules_dir).toSorted() : []
  // dev servers may skip the docs, but a build without them would ship an empty /api
  if (building && files.length === 0) {
    throw new Error(`No API docs in ${modules_dir}: run scripts/make_api_docs.py first`)
  }
  const modules = await Promise.all(
    files.map(async (file) => {
      const source = readFileSync(`${modules_dir}/${file}`, `utf8`)
      return assert_ok(await markdown.render(source, { filename: file }))
    }),
  )
  return { modules }
}
