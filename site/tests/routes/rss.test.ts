import { MODELS } from '$lib'
import { GET } from '$routes/rss.xml/+server'
import * as pkg from '$site/package.json'
import { describe, expect, it, vi } from 'vitest'

// Mock the MODELS array with a simplified version for testing
vi.mock(`$lib`, async () => {
  const actual = await vi.importActual(`$lib`)
  return {
    ...actual,
    MODELS: [
      {
        model_name: `Test Model`,
        model_key: `test-model`,
        model_params: 1000000,
        date_added: `2023-01-01`,
        date_published: `2023-01-15`,
        targets: `EFS_G`,
        model_type: `UIP`,
        hyperparams: {
          batch_size: 300,
          optimizer: `Adam`,
          loss: `MAE`,
          loss_weights: { energy: 10.0, force: 1.0 },
        },
        license: { code: `Apache-2.0` },
        training_set: [`MP 2022`],
        authors: [{ name: `Test Author`, affiliation: `Test University` }],
        metrics: {
          discovery: {
            full_test_set: {
              F1: 0.8,
              Precision: 0.85,
              Accuracy: 0.9,
              TPR: 0.7,
              TNR: 0.95,
              R2: 0.75,
              MAE: 0.1,
            },
          },
        },
        paper: `https://example.com/paper`,
        repo: `https://github.com/example/repo`,
      },
    ],
  }
})

// Mock the format_train_set function
vi.mock(`$lib/metrics`, () => ({
  format_train_set: () => `Test Training Set (1000 structures)`,
}))

describe(`RSS feed endpoint`, () => {
  it(`should return response with correct content type`, async () => {
    const response = await GET()
    expect(response.headers.get(`Content-Type`)).toBe(`application/xml`)
  })

  it(`should return valid XML with expected structure`, async () => {
    const response = await GET()
    const xml = await response.text()

    // Check RSS basics
    expect(xml).toContain(`<rss xmlns:atom="http://www.w3.org/2005/Atom" version="2.0">`)
    expect(xml).toContain(`<channel>`)
    expect(xml).toContain(`</channel>`)
    expect(xml).toContain(`</rss>`)

    // Check site info
    expect(xml).toContain(`<title>${pkg.name} - Model Leaderboard</title>`)
    expect(xml).toContain(`<description>${pkg.description}</description>`)

    // Check item structure
    expect(xml).toContain(`<item>`)
    expect(xml).toContain(`<title>Test Model</title>`)
    expect(xml).toContain(`<description><![CDATA[`)
    expect(xml).toContain(`</description>`)
    expect(xml).toContain(`<pubDate>`)
    expect(xml).toContain(`</item>`)
  })

  it(`should include model details in the description`, async () => {
    const response = await GET()
    const xml = await response.text()

    // Check for expected model details
    expect(xml).toContain(`<h2>Test Model</h2>`) // No "Model:" prefix
    expect(xml).toContain(`<strong>Authors:</strong> Test Author (Test University)`)
    expect(xml).toContain(`<strong>Date Published:</strong> 2023-01-15`)
    expect(xml).toContain(`<strong>Date Added:</strong> 2023-01-01`)
    expect(xml).toContain(`<strong>Training Set:</strong> Test Training Set`)
    expect(xml).toContain(`<strong>Parameters:</strong>`)
    expect(xml).toContain(`<strong>Model Type:</strong> UIP`)
    expect(xml).toContain(`<strong>Targets:</strong> EFS_G`)
    expect(xml).toContain(`<strong>Metrics:</strong>`)
    expect(xml).toContain(`F1: 0.8`)
    expect(xml).toContain(`R2: 0.75`)
    expect(xml).toContain(`<strong>Key Hyperparameters:</strong>`)
    expect(xml).toContain(`optimizer: Adam`)
    expect(xml).toContain(`loss: MAE`)
    expect(xml).toContain(`<strong>License:</strong> Apache-2.0`)
  })

  it(`should include links to model resources`, async () => {
    const response = await GET()
    const xml = await response.text()

    const base_url = pkg.homepage.endsWith(`/`) ? pkg.homepage : `${pkg.homepage}/`

    // Check for links
    expect(xml).toContain(`<link>${base_url}models/test-model</link>`)
    expect(xml).toContain(`<guid isPermaLink="true">${base_url}models/test-model</guid>`)
    expect(xml).toContain(`href="https://example.com/paper">Read paper</a>`)
    expect(xml).toContain(
      `href="https://github.com/example/repo">View code repository</a>`,
    )
  })

  it(`should sort models by date in descending order`, async () => {
    // Create a version of MODELS with multiple entries with different dates
    const multi_date_models = [
      {
        ...MODELS[0],
        model_name: `Newer Model`,
        model_key: `newer-model`,
        date_added: `2023-02-01`,
      },
      {
        ...MODELS[0],
        model_name: `Older Model`,
        model_key: `older-model`,
        date_added: `2022-12-01`,
      },
    ]

    // Mock MODELS with multiple entries
    const original_models = MODELS
    vi.mocked(MODELS).splice(0, MODELS.length, ...multi_date_models)

    const response = await GET()
    const xml = await response.text()

    // Newer model should appear before older model
    const newer_index = xml.indexOf(`Newer Model`)
    const older_index = xml.indexOf(`Older Model`)

    expect(newer_index).not.toBe(-1)
    expect(older_index).not.toBe(-1)
    expect(newer_index).toBeLessThan(older_index)

    // Restore original MODELS
    vi.mocked(MODELS).splice(0, MODELS.length, ...original_models)
  })
})
