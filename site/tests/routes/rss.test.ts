import { MODELS } from '$lib'
import { GET } from '$routes/rss.xml/+server'
import * as pkg from '$site/package.json'
import { describe, expect, it } from 'vitest'

describe(`RSS feed endpoint`, () => {
  it(`should return response with correct content type`, async () => {
    const response = await GET({ url: new URL(`https://example.com`) })
    expect(response.headers.get(`Content-Type`)).toBe(`application/xml`)
    expect(response.status).toBe(200)
  })

  it(`should return valid XML with expected structure`, async () => {
    const response = await GET({ url: new URL(`https://example.com`) })
    const xml = await response.text()

    // Check RSS basics with exact matches
    expect(xml).toMatch(
      /<rss xmlns:atom="http:\/\/www\.w3\.org\/2005\/Atom" version="2\.0">/,
    )
    expect(xml).toMatch(/<channel>[\s\S]*<\/channel>/)
    expect(xml).toContain(`</rss>`)

    // Check site info
    expect(xml).toContain(`<title>${pkg.name} - Model Leaderboard</title>`)
    expect(xml).toContain(`<description>${pkg.description}</description>`)

    // Check for required RSS elements
    expect(xml).toMatch(/<atom:link href=.*rel="self" type="application\/rss\+xml"\/>/)
    expect(xml).toMatch(/<link>.*<\/link>/)

    // Validate item structure
    const items = xml.match(/<item>[\s\S]*?<\/item>/g)
    expect(items).not.toBeNull()
    expect(items?.length).toBeGreaterThan(0)

    const first_item = items?.[0] || ``
    expect(first_item).toMatch(/<title>.*<\/title>/)
    expect(first_item).toMatch(/<link>.*<\/link>/)
    expect(first_item).toMatch(/<description><!\[CDATA\[[\s\S]*?\]\]><\/description>/)
    expect(first_item).toMatch(/<pubDate>.*<\/pubDate>/)
    expect(first_item).toMatch(/<guid isPermaLink="true">.*<\/guid>/)
  })

  it(`should include model details in the description`, async () => {
    // Skip test if no models available
    if (MODELS.length === 0) {
      console.warn(`Skipping test: No models available`)
      return
    }

    const response = await GET({ url: new URL(`https://example.com`) })
    const xml = await response.text()

    // Extract the CDATA content from the first item
    const cdata_match = xml.match(
      /<description><!\[CDATA\[([\s\S]*?)\]\]><\/description>/,
    )
    expect(cdata_match).not.toBeNull()

    const cdata_content = cdata_match?.[1] || ``

    // Check for expected model details with more specific patterns
    expect(cdata_content).toMatch(/<h2>[^<]+<\/h2>/) // Model name heading
    expect(cdata_content).toMatch(/<strong>Authors:<\/strong>[^<]+/)
    expect(cdata_content).toMatch(/<strong>Date Added:<\/strong>[^<]+/)
    expect(cdata_content).toMatch(/<strong>Training Set:<\/strong>[^<]+/)
    expect(cdata_content).toMatch(/<strong>Parameters:<\/strong>[^<]+/)
    expect(cdata_content).toMatch(/<strong>Model Type:<\/strong>[^<]+/)

    // Check that required sections exist
    expect(cdata_content).toMatch(/<strong>Metrics:<\/strong>/)

    // Check for proper HTML structure
    const strongTags = cdata_content.match(/<strong>/g)
    const strongCloseTags = cdata_content.match(/<\/strong>/g)
    const divTags = cdata_content.match(/<div>/g)
    const divCloseTags = cdata_content.match(/<\/div>/g)

    expect(strongTags?.length ?? 0).toBeGreaterThanOrEqual(5)
    expect(strongCloseTags?.length ?? 0).toBeGreaterThanOrEqual(5)
    expect(divTags?.length ?? 0).toBeGreaterThanOrEqual(0)
    expect(divCloseTags?.length ?? 0).toBeGreaterThanOrEqual(0)
  })

  it(`should include links to model resources`, async () => {
    // Skip test if no models available
    if (MODELS.length === 0) {
      console.warn(`Skipping test: No models available`)
      return
    }

    const response = await GET({ url: new URL(`https://example.com`) })
    const xml = await response.text()

    const base_url = pkg.homepage.endsWith(`/`) ? pkg.homepage : `${pkg.homepage}/`

    // Check for model links with exact pattern matching
    expect(xml).toMatch(new RegExp(`<link>${base_url}models/[^<]+</link>`))
    expect(xml).toMatch(
      new RegExp(`<guid isPermaLink="true">${base_url}models/[^<]+</guid>`),
    )

    // Check for at least one model from the real data
    if (MODELS.length > 0) {
      const model = MODELS[0]
      expect(xml).toContain(model.model_name)
      expect(xml).toContain(`models/${model.model_key}`)
    }

    // Check for links to paper and repo in the description if available
    const cdata_match = xml.match(
      /<description><!\[CDATA\[([\s\S]*?)\]\]><\/description>/,
    )
    const cdata_content = cdata_match?.[1] || ``

    // Check for link patterns
    if (cdata_content.includes(`paper`)) {
      expect(cdata_content).toMatch(/href="https?:\/\/[^"]+">[^<]+<\/a>/)
    }
  })

  it(`should sort models by date in descending order`, async () => {
    // Skip test if not enough models available
    if (MODELS.length < 2) {
      console.warn(`Skipping test: Not enough models available to test sorting`)
      return
    }

    const response = await GET({ url: new URL(`https://example.com`) })
    const xml = await response.text()

    // Find models with different dates
    const sorted_models = [...MODELS].sort(
      (a, b) => new Date(b.date_added).getTime() - new Date(a.date_added).getTime(),
    )

    // Take the first two models with different dates
    let newer_model = null
    let older_model = null

    // Use for...of loop instead of index-based loop
    let prev_model = sorted_models[0]
    for (const current_model of sorted_models.slice(1)) {
      if (prev_model.date_added !== current_model.date_added) {
        newer_model = prev_model
        older_model = current_model
        break
      }
      prev_model = current_model
    }

    // Skip test if we couldn't find two models with different dates
    if (!newer_model || !older_model) {
      console.warn(`Skipping test: Couldn't find two models with different dates`)
      return
    }

    // Verify the models have different dates
    expect(newer_model.date_added).not.toBe(older_model.date_added)
    expect(new Date(newer_model.date_added).getTime()).toBeGreaterThan(
      new Date(older_model.date_added).getTime(),
    )

    // Newer model should appear before older model in the XML
    const newer_index = xml.indexOf(newer_model.model_name)
    const older_index = xml.indexOf(older_model.model_name)

    expect(newer_index).not.toBe(-1)
    expect(older_index).not.toBe(-1)
    expect(newer_index).toBeLessThan(older_index)

    // If test fails, this would help with debugging
    if (newer_index === -1) {
      console.warn(`Newer model '${newer_model.model_name}' not found in RSS feed`)
    }
    if (older_index === -1) {
      console.warn(`Older model '${older_model.model_name}' not found in RSS feed`)
    }
    if (newer_index >= older_index) {
      console.warn(
        `Newer model (${newer_model.model_name}, ${newer_model.date_added}) should appear before older model (${older_model.model_name}, ${older_model.date_added})`,
      )
    }

    // Additional check: all items should have a pubDate in correct format
    const pub_dates = xml.match(/<pubDate>[^<]+<\/pubDate>/g) || []
    expect(pub_dates.length).toBeGreaterThan(0)
    for (const date_str of pub_dates) {
      const date_content = date_str.replace(/<\/?pubDate>/g, ``)
      // Check that this parses as a valid date
      expect(new Date(date_content).toString()).not.toBe(`Invalid Date`)
    }
  })

  it(`should have correct self-reference URL and absolute URLs in descriptions`, async () => {
    const mock_url = new URL(`https://example.com/rss.xml`)

    const response = await GET({ url: mock_url })
    const xml = await response.text()

    // Check that the self-reference URL matches document URL
    expect(xml).toContain(
      `<atom:link href="${pkg.homepage}/rss.xml" rel="self" type="application/rss+xml"/>`,
    )

    // Extract descriptions to check for relative URLs
    const descriptions =
      xml.match(/<description><!\[CDATA\[([\s\S]*?)\]\]><\/description>/g) || []
    expect(descriptions.length).toBeGreaterThan(0)

    // URLs in descriptions should be absolute
    const url_matches = (descriptions[0] || ``).match(/href="([^"]+)"/g) || []
    expect(url_matches.length).toBeGreaterThan(0)

    // All URLs should be absolute (start with https://)
    for (const url_match of url_matches) {
      const url = url_match.replace(/href="|"/g, ``)
      expect(url.startsWith(`https://`)).toBe(true)
    }
  })
})
