import { describe, expect, it, vi } from 'vitest'
import { entries, load } from '$routes/models/[slug]/+page.server'

vi.mock(`$lib`, () => ({
  MODELS: [{ model_key: `model-key` }],
}))
vi.mock(`$lib/server/md`, () => ({ read_md_per_system: vi.fn() }))

describe(`model-key routes`, () => {
  it(`resolves the canonical key`, async () => {
    expect(entries()).toEqual([{ slug: `model-key` }])
    const result = await load({
      params: { slug: `model-key` },
    } as Parameters<typeof load>[0])
    expect(result).toMatchObject({ model: { model_key: `model-key` } })
  })

  it(`rejects unknown keys`, async () => {
    await expect(
      load({
        params: { slug: `missing` },
      } as Parameters<typeof load>[0]),
    ).rejects.toMatchObject({ status: 404 })
  })
})
