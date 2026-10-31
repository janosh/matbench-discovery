import { default_highlighter } from 'svelte-widgets/highlight'
import { assert_ok, create_markdown } from 'svelte-widgets/markdown'

// Generated from docstrings by scripts/make_api_docs.py. Rendering as plain Markdown
// instead of compiling Svelte keeps docstring braces like `({} if absent)` literal.
const sources = import.meta.glob<string>(`./modules/*.md`, {
  query: `?raw`,
  import: `default`,
  eager: true,
})

export const load = async () => {
  const markdown = create_markdown({ highlight: default_highlighter.highlight })
  const entries = Object.entries(sources).toSorted(([path_a], [path_b]) =>
    path_a.localeCompare(path_b),
  )
  const modules = await Promise.all(
    entries.map(async ([filename, source]) =>
      assert_ok(await markdown.render(source, { filename })),
    ),
  )
  return { modules }
}
