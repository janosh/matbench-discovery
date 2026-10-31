import { by_benchmark_added_desc } from '#lib'
import { MODELS } from '#lib/models.svelte.js'
import { GET } from '#routes/rss.xml/+server.js'
import pkg from '#site/package.json'
import { format_num } from 'matterviz/labels'
import { escape_html } from 'svelte-widgets/utils'
import { describe, expect, it, onTestFinished } from 'vite-plus/test'

const base_url = pkg.homepage.endsWith(`/`) ? pkg.homepage : `${pkg.homepage}/`
const sorted_models = MODELS.toSorted(by_benchmark_added_desc)
const cdata_descriptions = (xml: string): string[] =>
  [
    ...xml.matchAll(/<description><!\[CDATA\[(?<cdata>[\s\S]*?)\]\]><\/description>/g),
  ].map((match) => match.groups?.cdata ?? ``)

describe(`RSS feed endpoint`, () => {
  it(`returns an RSS document with one linked item per model, newest first`, async () => {
    const response = GET()
    expect(response.headers.get(`Content-Type`)).toBe(`application/xml`)
    expect(response.status).toBe(200)
    const xml = await response.text()

    expect(xml).toMatch(
      /^<rss xmlns:atom="http:\/\/www\.w3\.org\/2005\/Atom" version="2\.0">\s*<channel>[\s\S]*<\/channel>\s*<\/rss>$/,
    )
    expect(xml).toContain(`<title>${escape_html(pkg.name)}</title>`)
    expect(xml).toContain(`<description>${escape_html(pkg.description)}</description>`)
    expect(xml).toContain(
      `<atom:link href="${base_url}rss.xml" rel="self" type="application/rss+xml"/>`,
    )

    const items = [...xml.matchAll(/<item>[\s\S]*?<\/item>/g)].map((match) => match[0])
    expect(items).toHaveLength(MODELS.length)
    for (const [idx, item] of items.entries()) {
      const { model_key, model_name } = sorted_models[idx]
      const model_url = `${base_url}models/${model_key}`
      expect(/<title>(?<title>[^<]*)<\/title>/.exec(item)?.groups?.title).toBe(
        escape_html(model_name),
      )
      expect(item).toContain(`<link>${model_url}</link>`)
      expect(item).toContain(`<guid isPermaLink="true">${model_url}</guid>`)
      const pub_date = /<pubDate>(?<date>[^<]+)<\/pubDate>/.exec(item)?.groups?.date
      expect(new Date(pub_date ?? ``).toString()).not.toBe(`Invalid Date`)
    }
  })

  it(`escapes XML special characters in item titles`, async () => {
    const [model] = MODELS
    const { model_name } = model
    model.model_name = `A & <B>`
    onTestFinished(() => {
      model.model_name = model_name
    })
    const xml = await GET().text()
    expect(xml).toContain(`<title>A &amp; &lt;B&gt;</title>`)
    expect(xml).not.toContain(`<title>A & <B></title>`)
  })

  it(`describes every model with ordered details and absolute resource links`, async () => {
    const descriptions = cdata_descriptions(await GET().text())
    expect(descriptions).toHaveLength(sorted_models.length)

    const [first] = descriptions
    expect(first).toMatch(/<h2>[^<]+<\/h2>/)
    expect(first).toMatch(/<strong>Parameters:<\/strong>[^<]+/)
    expect(first).toContain(`<strong>Role:</strong>`)
    expect(first.indexOf(`<strong>Metrics:</strong>`)).toBeGreaterThan(-1)
    expect(first.indexOf(`<strong>Metrics:</strong>`)).toBeLessThan(
      first.indexOf(`<strong>Authors:</strong>`),
    )

    // each description links to its detail page plus paper and repo when set (some
    // models, e.g. EMA-GNN, have paper: null); all hrefs are absolute
    expect(sorted_models.some((model_data) => model_data.paper)).toBe(true)
    for (const [idx, model_data] of sorted_models.entries()) {
      // metrics quote the leaderboard's (CPS) discovery set, not the full test set
      const { unique_prototypes, full_test_set } = model_data.metrics?.discovery ?? {}
      const f1 = (unique_prototypes ?? full_test_set)?.F1
      if (f1 != null) expect(descriptions[idx]).toContain(`F1: ${format_num(f1)},`)
      const links = [...descriptions[idx].matchAll(/<a href="[^"]+">[^<]+<\/a>/g)].map(
        (match) => match[0],
      )
      expect(links, model_data.model_key).toStrictEqual([
        `<a href="${base_url}models/${model_data.model_key}">View model details</a>`,
        ...(model_data.paper ? [`<a href="${model_data.paper}">Read paper</a>`] : []),
        ...(model_data.repo
          ? [`<a href="${model_data.repo}">View code repository</a>`]
          : []),
      ])
      for (const match of descriptions[idx].matchAll(/href="(?<url>[^"]+)"/g))
        expect(match.groups?.url).toMatch(/^https:\/\//)
    }
  })

  it(`should include underscore-containing hyperparams in descriptions`, async () => {
    // regression: the filter used !key.includes('_') instead of startsWith, dropping
    // nearly every hyperparam (max_steps, ase_optimizer, ...) so the section was empty
    const inner_underscore = (key: string) => key.includes(`_`) && !key.startsWith(`_`)
    const hyperparameter_key = MODELS.flatMap((model) =>
      Object.entries(model.hyperparams ?? {}).flatMap(([namespace, values]) =>
        Object.keys(values).map((key) => `${namespace}.${key}`),
      ),
    ).find(inner_underscore)
    if (!hyperparameter_key)
      throw new Error(`no model with underscore-containing hyperparams`)

    const xml = await GET().text()
    expect(xml).toContain(`${hyperparameter_key}: `)
  })
})
