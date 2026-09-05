import { by_benchmark_added_desc, MODELS } from '$lib'
import { GET } from '$routes/rss.xml/+server'
import pkg from '$site/package.json'
import { describe, expect, it } from 'vitest'

describe(`RSS feed endpoint`, () => {
  const extract_first_cdata = (xml: string): string => {
    const match = /<description><!\[CDATA\[(?<cdata>[\s\S]*?)\]\]><\/description>/.exec(
      xml,
    )
    return match?.groups?.cdata ?? ``
  }

  it(`should return valid XML with expected structure and content type`, async () => {
    const response = GET()
    expect(response.headers.get(`Content-Type`)).toBe(`application/xml`)
    expect(response.status).toBe(200)
    const xml = await response.text()

    expect(xml).toMatch(
      /<rss xmlns:atom="http:\/\/www\.w3\.org\/2005\/Atom" version="2\.0">/,
    )
    expect(xml).toMatch(/<channel>[\s\S]*<\/channel>/)
    expect(xml).toContain(`</rss>`)

    expect(xml).toContain(`<title>${pkg.name}</title>`)
    expect(xml).toContain(`<description>${pkg.description}</description>`)

    expect(xml).toMatch(/<atom:link href=.*rel="self" type="application\/rss\+xml"\/>/)
    expect(xml).toMatch(/<link>.*<\/link>/)

    // Validate item structure: the feed contains one item per model
    const items = [...xml.matchAll(/<item>[\s\S]*?<\/item>/g)]
    expect(items).toHaveLength(MODELS.length)

    const first_item = items[0]?.[0] ?? ``
    expect(first_item).toMatch(/<title>.*<\/title>/)
    expect(first_item).toMatch(/<link>.*<\/link>/)
    expect(first_item).toMatch(/<description><!\[CDATA\[[\s\S]*?\]\]><\/description>/)
    expect(first_item).toMatch(/<pubDate>.*<\/pubDate>/)
    expect(first_item).toMatch(/<guid isPermaLink="true">.*<\/guid>/)
  })

  it(`should include model details in the description`, async () => {
    const response = GET()
    const xml = await response.text()

    const cdata_content = extract_first_cdata(xml)
    expect(cdata_content).not.toBe(``)

    expect(cdata_content).toMatch(/<h2>[^<]+<\/h2>/) // Model name heading
    expect(cdata_content).toMatch(/<strong>Metrics:<\/strong>/)
    expect(cdata_content).toMatch(/<strong>Parameters:<\/strong>[^<]+/)

    // Role is present for every model.
    expect(cdata_content).toContain(`<strong>Role:</strong>`)

    // Authors section should appear later in the content
    const metrics_pos = cdata_content.indexOf(`<strong>Metrics:</strong>`)
    const authors_pos = cdata_content.indexOf(`<strong>Authors:</strong>`)
    expect(metrics_pos).toBeLessThan(authors_pos)

    const strong_tags = [...cdata_content.matchAll(/<strong>/g)]
    const strong_close_tags = [...cdata_content.matchAll(/<\/strong>/g)]

    expect(strong_tags.length).toBeGreaterThanOrEqual(5)
    expect(strong_close_tags.length).toBeGreaterThanOrEqual(5)
  })

  it(`should include links to model resources`, async () => {
    const response = GET()
    const xml = await response.text()

    const base_url = pkg.homepage.endsWith(`/`) ? pkg.homepage : `${pkg.homepage}/`

    expect(xml).toMatch(new RegExp(`<link>${base_url}models/[^<]+</link>`))
    expect(xml).toMatch(
      new RegExp(`<guid isPermaLink="true">${base_url}models/[^<]+</guid>`),
    )

    const model = MODELS[0]
    expect(xml).toContain(model.model_name)
    expect(xml).toContain(`models/${model.model_key}`)

    // every item's description links exactly to its detail page plus the model's
    // paper and repo when set (some models, e.g. EMA-GNN, have paper: null)
    const descriptions = [
      ...xml.matchAll(/<description><!\[CDATA\[(?<cdata>[\s\S]*?)\]\]><\/description>/g),
    ].map((match) => match.groups?.cdata ?? ``)
    const sorted_models = MODELS.toSorted(by_benchmark_added_desc)
    expect(descriptions).toHaveLength(sorted_models.length)
    expect(sorted_models.some((model_data) => model_data.paper)).toBe(true)

    for (const [idx, model_data] of sorted_models.entries()) {
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

  it(`should sort models by date in descending order`, async () => {
    const xml = await GET().text()

    const expected_names = MODELS.toSorted(by_benchmark_added_desc).map(
      (model) => model.model_name,
    )
    const item_names = [...xml.matchAll(/<item>\s*<title>(?<name>[^<]+)<\/title>/g)].map(
      (match) => match.groups?.name ?? ``,
    )
    expect(item_names).toStrictEqual(expected_names)

    // every item carries a parseable pubDate
    const pub_dates = [...xml.matchAll(/<pubDate>(?<date>[^<]+)<\/pubDate>/g)]
    expect(pub_dates).toHaveLength(MODELS.length)
    for (const match of pub_dates) {
      expect(new Date(match.groups?.date ?? ``).toString()).not.toBe(`Invalid Date`)
    }
  })

  it(`should have correct self-reference URL and absolute URLs in descriptions`, async () => {
    const response = GET()
    const xml = await response.text()

    expect(xml).toContain(
      `<atom:link href="${pkg.homepage}/rss.xml" rel="self" type="application/rss+xml"/>`,
    )

    const descriptions = [
      ...xml.matchAll(/<description><!\[CDATA\[(?<cdata>[\s\S]*?)\]\]><\/description>/g),
    ]
    expect(descriptions).toHaveLength(MODELS.length)

    // URLs in descriptions should be absolute
    const description = descriptions[0]?.groups?.cdata ?? ``
    const url_matches = [...description.matchAll(/href="(?<url>[^"]+)"/g)]
    expect(url_matches.length).toBeGreaterThan(0)

    // All URLs should be absolute (start with https://)
    for (const match of url_matches) {
      expect(match.groups?.url?.startsWith(`https://`)).toBe(true)
    }
  })
})
