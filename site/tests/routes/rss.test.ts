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

  it(`should return response with correct content type`, () => {
    const response = GET()
    expect(response.headers.get(`Content-Type`)).toBe(`application/xml`)
    expect(response.status).toBe(200)
  })

  it(`should return valid XML with expected structure`, async () => {
    const response = GET()
    const xml = await response.text()

    // Check RSS basics with exact matches
    expect(xml).toMatch(
      /<rss xmlns:atom="http:\/\/www\.w3\.org\/2005\/Atom" version="2\.0">/,
    )
    expect(xml).toMatch(/<channel>[\s\S]*<\/channel>/)
    expect(xml).toContain(`</rss>`)

    // Check site info
    expect(xml).toContain(`<title>${pkg.name}</title>`)
    expect(xml).toContain(`<description>${pkg.description}</description>`)

    // Check for required RSS elements
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

    // Extract the CDATA content from the first item
    const cdata_content = extract_first_cdata(xml)
    expect(cdata_content).not.toBe(``)

    // Check for expected model details in the correct order
    expect(cdata_content).toMatch(/<h2>[^<]+<\/h2>/) // Model name heading
    expect(cdata_content).toMatch(/<strong>Metrics:<\/strong>/)
    expect(cdata_content).toMatch(/<strong>Parameters:<\/strong>[^<]+/)

    // Role is present for every model.
    expect(cdata_content).toContain(`<strong>Role:</strong>`)

    // Authors section should appear later in the content
    const metrics_pos = cdata_content.indexOf(`<strong>Metrics:</strong>`)
    const authors_pos = cdata_content.indexOf(`<strong>Authors:</strong>`)
    expect(metrics_pos).toBeLessThan(authors_pos)

    // Check for proper HTML structure
    const strong_tags = [...cdata_content.matchAll(/<strong>/g)]
    const strong_close_tags = [...cdata_content.matchAll(/<\/strong>/g)]

    expect(strong_tags.length).toBeGreaterThanOrEqual(5)
    expect(strong_close_tags.length).toBeGreaterThanOrEqual(5)
  })

  it(`should include links to model resources`, async () => {
    const response = GET()
    const xml = await response.text()

    const base_url = pkg.homepage.endsWith(`/`) ? pkg.homepage : `${pkg.homepage}/`

    // Check for model links with exact pattern matching
    expect(xml).toMatch(new RegExp(`<link>${base_url}models/[^<]+</link>`))
    expect(xml).toMatch(
      new RegExp(`<guid isPermaLink="true">${base_url}models/[^<]+</guid>`),
    )

    // Check for at least one model from the real data
    const model = MODELS[0]
    expect(xml).toContain(model.model_name)
    expect(xml).toContain(`models/${model.model_key}`)

    // Check for links to paper and repo in the description if available
    const cdata_content = extract_first_cdata(xml)

    // Check for link patterns
    expect(
      !cdata_content.includes(`paper`) ||
        /href="https?:\/\/[^"]+">[^<]+<\/a>/.test(cdata_content),
    ).toBe(true)
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

    // Check that the self-reference URL matches document URL
    expect(xml).toContain(
      `<atom:link href="${pkg.homepage}/rss.xml" rel="self" type="application/rss+xml"/>`,
    )

    // Extract descriptions to check for relative URLs
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
