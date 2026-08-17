import { describe, expect, it, vi } from 'vitest'
import { entries, load } from '$routes/models/[slug]/+page.server'

vi.mock(`$lib`, () => ({
  MODELS: [{ model_key: `new-key`, model_key_aliases: [`old-key`] }],
}))
vi.mock(`$lib/server/md`, () => ({ read_md_per_system: vi.fn() }))

describe(`model-key aliases`, () => {
  it(`generates alias routes that permanently redirect to the canonical key`, async () => {
    expect(entries()).toEqual([{ slug: `new-key` }, { slug: `old-key` }])
    await expect(
      load({
        params: { slug: `old-key` },
        url: new URL(`https://example.com/models/old-key?energy_tab=each`),
      } as Parameters<typeof load>[0]),
    ).rejects.toMatchObject({
      status: 308,
      location: `/models/new-key?energy_tab=each`,
    })
  })

  it(`resolves the canonical key`, async () => {
    const result = await load({
      params: { slug: `new-key` },
      url: new URL(`https://example.com/models/new-key`),
    } as Parameters<typeof load>[0])
    expect(result).toMatchObject({ model: { model_key: `new-key` } })
  })

  it(`rejects unknown keys`, async () => {
    await expect(
      load({
        params: { slug: `missing` },
        url: new URL(`https://example.com/models/missing`),
      } as Parameters<typeof load>[0]),
    ).rejects.toMatchObject({ status: 404 })
  })
})
