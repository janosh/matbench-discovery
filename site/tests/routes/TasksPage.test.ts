import MODELING_TASKS from '$pkg/modeling-tasks.yml'
import TasksPage from '$routes/tasks/+page.svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

describe(`Tasks Page`, () => {
  it(`renders a subpage card per task page linking to its route`, () => {
    mount(TasksPage, { target: document.body })

    expect(doc_query(`.subpage-grid h1`).textContent).toBe(`Benchmark Tasks`)
    const cards = [
      ...document.querySelectorAll<HTMLAnchorElement>(`.subpage-grid a.card`),
    ]
    // CPS has no dedicated page
    const expected = [`discovery`, `geo_opt`, `phonons`, `diatomics`, `md`] as const
    expect(cards.map((card) => card.getAttribute(`href`))).toEqual(
      expected.map((key) => `/tasks/${key.replaceAll(`_`, `-`)}`),
    )
    expect(cards.map((card) => card.querySelector(`h2`)?.textContent)).toEqual(
      expected.map((key) => MODELING_TASKS[key].label),
    )
    for (const [idx, key] of expected.entries()) {
      expect(cards[idx].textContent).toContain(MODELING_TASKS[key].description)
      expect(cards[idx].querySelector(`svg`)).not.toBeNull()
    }
  })
})
