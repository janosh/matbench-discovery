import MODELING_TASKS from '$pkg/modeling-tasks.yml'
import BenchmarksPage from '$routes/benchmarks/+page.svelte'
import { handle } from '$site/src/hooks.server'
import { describe, expect, it, vi } from 'vitest'
import { doc_query, mount } from '../index'

describe(`Benchmarks Page`, () => {
  it.each([``, `?cps_weights=1,0,0`])(
    `redirects /tasks%s in the server handler with HTTP 307`,
    async (suffix) => {
      const resolve = vi.fn()
      await expect(
        handle({
          event: {
            url: new URL(`http://localhost/tasks${suffix}`),
          } as Parameters<typeof handle>[0][`event`],
          resolve,
        }),
      ).rejects.toMatchObject({
        status: 307,
        location: `/benchmarks${suffix}`,
      })
      expect(resolve).not.toHaveBeenCalled()
    },
  )

  it(`merges tasks and datasets with distinct links to results and test sets`, () => {
    mount(BenchmarksPage, { target: document.body })

    expect(doc_query(`h1`).textContent).toBe(`Benchmarks`)
    const cards = [
      ...document.querySelectorAll<HTMLAnchorElement>(`.benchmark-grid h2 a`),
    ]
    // CPS has no dedicated page
    const expected = [`discovery`, `geo_opt`, `phonons`, `md`, `diatomics`] as const
    expect(cards.map((card) => card.getAttribute(`href`))).toEqual(
      expected.map((key) => `/benchmarks/${key.replaceAll(`_`, `-`)}`),
    )
    expect(cards.map((card) => card.textContent?.trim())).toEqual(
      expected.map((key) => MODELING_TASKS[key].label),
    )
    expect(document.querySelectorAll(`.benchmark-grid .periodic-table`)).toHaveLength(5)
    expect(document.querySelector(`.data-files-list`)).not.toBeNull()
    expect(doc_query(`a[href="/data/sets"]`)).not.toBeNull()
    for (const [idx, key] of expected.entries()) {
      const card = doc_query(
        `article:nth-child(${idx + 1})`,
        doc_query(`.benchmark-grid`),
      )
      const table = doc_query(`h2 + figure .periodic-table`, card)
      const tiles = [...table.querySelectorAll<HTMLElement>(`[data-element-symbol]`)]
      expect(tiles.filter((tile) => tile.style.opacity !== `0.15`)).toHaveLength(
        [85, 85, 34, 22, 92][idx],
      )
      expect(table.querySelector(`.element-tile:not([data-element-symbol])`)).toBeNull()
      expect(doc_query(`.colorbar`, table).textContent?.includes(`(log)`)).toBe(
        key !== `diatomics`,
      )
      expect(card.textContent).toContain(MODELING_TASKS[key].description)
      expect(cards[idx].querySelector(`svg`)).not.toBeNull()
      expect(
        card.querySelector(`a[href="/benchmarks/${key.replaceAll(`_`, `-`)}#test-set"]`),
      ).not.toBeNull()
    }
  })
})
