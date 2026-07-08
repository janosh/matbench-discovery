import { type Label, TableControls } from '$lib'
import { ALL_TRAINING_SETS, make_table_filters } from '$lib/models.svelte'
import { OPENNESS_OPTIONS } from '$lib/url-state.svelte'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { doc_query, mount } from '../index'

describe(`TableControls`, () => {
  const sample_columns: Label[] = [
    { key: `model`, label: `Model`, description: `Model name`, visible: true },
    { key: `f1`, label: `F1`, description: `F1 Score`, visible: true },
    { key: `daf`, label: `DAF`, description: `DAF Score`, visible: true },
    { key: `rmse`, label: `RMSE`, description: `RMSE`, visible: false },
  ]

  const summary_for = (text: string): HTMLElement => {
    const summary = [...document.querySelectorAll(`details.filter-menu summary`)].find(
      (el) => el.textContent?.includes(text),
    )
    if (!summary) throw new Error(`No filter summary found containing: ${text}`)
    return summary as HTMLElement
  }

  it(`renders filter dropdowns and heatmap toggle`, () => {
    mount(TableControls, { target: document.body })

    expect(summary_for(`Training data`)).toBeDefined()
    expect(summary_for(`Openness`)).toBeDefined()
    const labels = [...document.querySelectorAll(`label`)].map((label) =>
      label.textContent?.replaceAll(/\s+/g, ` `).trim(),
    )
    expect(labels).toContain(`Heatmap`)
  })

  it(`training-data dropdown lists all datasets with only/not checkboxes`, async () => {
    const filters = make_table_filters()
    mount(TableControls, { target: document.body, props: { filters } })
    await tick()

    const dropdown = summary_for(`Training data`).closest(`details`)
    const rows = dropdown?.querySelectorAll(`.filter-row`) ?? []
    expect(rows).toHaveLength(ALL_TRAINING_SETS.length)

    // check `only` for the first dataset: require-mode filter becomes active,
    // its checkbox checks, and the summary shows a count badge
    const first_row = rows[0]
    const [only_box, not_box] = first_row.querySelectorAll<HTMLInputElement>(`input`)
    only_box.click()
    await tick()
    expect(filters.training[ALL_TRAINING_SETS[0]]).toBe(`require`)
    expect(only_box.checked).toBe(true)
    expect(summary_for(`Training data (1)`)).toBeDefined()

    // checking `not` on the same dataset flips the mode (mutually exclusive)
    not_box.click()
    await tick()
    expect(filters.training[ALL_TRAINING_SETS[0]]).toBe(`exclude`)
    expect(only_box.checked).toBe(false)
    expect(not_box.checked).toBe(true)

    // clear-filters button resets everything
    doc_query<HTMLButtonElement>(`button.clear-filters`).click()
    await tick()
    expect(filters.n_active).toBe(0)
  })

  it(`openness dropdown toggles values but never hides the last one`, async () => {
    const filters = make_table_filters()
    mount(TableControls, { target: document.body, props: { filters } })
    await tick()

    const dropdown = summary_for(`Openness`).closest(`details`)
    const boxes = dropdown?.querySelectorAll<HTMLInputElement>(`input`) ?? []
    expect(boxes).toHaveLength(OPENNESS_OPTIONS.length)
    expect([...boxes].every((box) => box.checked)).toBe(true)

    boxes[1].click() // hide OSCD
    await tick()
    expect(filters.openness).toStrictEqual([`OSOD`, `CSOD`, `CSCD`])
    expect(summary_for(`Openness (3/4)`)).toBeDefined()
  })

  it(`opens and closes column visibility panel`, async () => {
    mount(TableControls, {
      target: document.body,
      props: { columns: sample_columns },
    })
    await tick()

    const toggle_btn = doc_query(`.column-toggles summary`)
    const details = doc_query<HTMLDetailsElement>(`.column-toggles`)
    expect(toggle_btn.textContent?.trim()).toBe(`Columns`)

    // Should start closed
    expect(details.open).toBe(false)

    // Open panel
    toggle_btn.click()
    expect(details.open).toBe(true)

    // Verify column menu is visible
    expect(
      doc_query(`.column-menu`).querySelectorAll(`input[type="checkbox"]`),
    ).toHaveLength(sample_columns.length)

    // Close panel
    toggle_btn.click()
    expect(details.open).toBe(false)
  })

  it(`toggles column visibility checkboxes`, async () => {
    mount(TableControls, {
      target: document.body,
      props: { columns: [...sample_columns] },
    })
    await tick()

    const toggle_btn = doc_query(`.column-toggles summary`)
    toggle_btn.click()

    const column_menu = doc_query(`.column-menu`)
    expect(column_menu.getAttribute(`role`)).toBe(`group`)

    const column_checkboxes =
      column_menu.querySelectorAll<HTMLInputElement>(`input[type="checkbox"]`)
    expect(column_checkboxes).toHaveLength(sample_columns.length)
    const checkbox_labels = [...column_menu.querySelectorAll(`label`)].map((label) =>
      label.textContent?.trim(),
    )
    expect(checkbox_labels).toStrictEqual(sample_columns.map((column) => column.label))

    // Test toggling first checkbox (Model column, initially visible)
    const [first_checkbox] = column_checkboxes
    expect(first_checkbox.checked).toBe(true)

    first_checkbox.click()
    expect(first_checkbox.checked).toBe(false)

    first_checkbox.click()
    // Verify checkbox state changed
    expect(first_checkbox.checked).toBe(true)
  })
})
