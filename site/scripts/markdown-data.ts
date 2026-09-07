import { assert_ok, create_markdown, render_markdown } from 'svelte-widgets/markdown'

const engine = create_markdown({ raw_html: `omit` })

const record = (value: unknown, context: string): Record<string, unknown> => {
  if (!value || typeof value !== `object` || Array.isArray(value))
    throw new TypeError(`${context}: expected a mapping`)
  return value as Record<string, unknown>
}

async function render_text(value: unknown, context: string, suffix = ``) {
  if (typeof value !== `string`)
    throw new TypeError(`${context}: expected Markdown text, received ${typeof value}`)
  const document = assert_ok(
    await engine.parse(value + suffix, { filename: context, dialect: `markdown` }),
  )
  return assert_ok(await render_markdown(document))
}

// Enrich YAML data during Vite transforms, keeping parsing and mutation out of client modules.
export async function render_data_markdown(
  data: unknown,
  filename: string,
): Promise<unknown> {
  const file = filename.replaceAll(`\\`, `/`)
  const datasets = file.endsWith(`/data/datasets.yml`)
  const data_files = file.endsWith(`/matbench_discovery/data-files.yml`)
  const model = /\/models\/[^_][^/]*\/[^_][^/]*\.yml$/u.test(file)
  if (!datasets && !data_files && !model) return data
  const entries = record(data, filename)

  if (model) {
    if (entries.notes === undefined) return data
    const notes = record(entries.notes, `${filename}: notes`)
    const html = record((notes.html ??= {}), `${filename}: notes.html`)
    for (const [key, note] of Object.entries(notes)) {
      if (typeof note !== `string` || key in html) continue
      const rendered = await render_text(note, `${filename}: notes.${key}`)
      if (rendered) html[key] = rendered
    }
    return data
  }

  let references = ``
  if (data_files) {
    const { _links: links } = entries
    if (typeof links !== `string`)
      throw new TypeError(`${filename}: _links must be a string`)
    references = `\n\n${links}`
  }
  for (const [key, value] of Object.entries(entries)) {
    if (data_files && key.startsWith(`_`)) continue
    const entry = record(value, `${filename}: ${key}`)
    entry[datasets ? `description_html` : `html`] = await render_text(
      entry.description,
      `${filename}: ${key}.description`,
      references,
    )
    if (datasets) entry.slug = key.toLowerCase().replaceAll(/[\s_]+/g, `-`)
  }
  return data
}
