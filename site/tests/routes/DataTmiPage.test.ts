import DataTmiPage from '$routes/data/tmi/+page.svelte'
import { tick } from 'svelte'
import { describe, expect, it } from 'vitest'
import { checkbox_for, doc_query, mount_with_url } from '../index'

const checked_radio = (): HTMLInputElement | null =>
  document.querySelector(`input[type="radio"][name="filter"]:checked`)

const radio_values = (): string[] =>
  [
    ...document.querySelectorAll<HTMLInputElement>(`input[type="radio"][name="filter"]`),
  ].map((radio) => radio.value)

describe(`Data TMI Page`, () => {
  it(`defaults filter to first arity key and toggles unchecked`, async () => {
    await mount_with_url(DataTmiPage, `http://localhost/data/tmi`)

    const default_filter = radio_values().find((value) => value.startsWith(`arity=`))
    expect(checked_radio()?.value).toBe(default_filter)
    expect(checkbox_for(`Normalize by data set size`).checked).toBe(false)
    expect(checkbox_for(`Log color scale`).checked).toBe(false)
    // shared PtableHeatmap renders the count color bar with the filter in its title
    expect(document.querySelector(`.periodic-table .colorbar`)?.textContent).toContain(
      `WBM element counts for ${default_filter}`,
    )
  })

  it(`restores filter and toggles from URL params`, async () => {
    const target_filter = `batch=5`
    await mount_with_url(
      DataTmiPage,
      `http://localhost/data/tmi?filter=${encodeURIComponent(
        target_filter,
      )}&normalized=1&log=1&color_scale=interpolatePlasma`,
    )

    expect(checked_radio()?.value).toBe(target_filter)
    expect(checkbox_for(`Normalize by data set size`).checked).toBe(true)
    expect(checkbox_for(`Log color scale`).checked).toBe(true)
    expect(new URL(location.href).searchParams.get(`color_scale`)).toBe(
      `interpolatePlasma`,
    )
    const scale_input = doc_query<HTMLInputElement>(`input[aria-label="Color scale"]`)
    const picker = scale_input.closest(`.multiselect`)
    expect(picker?.querySelector(`ul.selected`)?.textContent).toContain(`Plasma`)
    scale_input.focus()
    await tick()
    const viridis_option = [
      ...(picker?.querySelectorAll<HTMLElement>(`ul.options li[aria-posinset]`) ?? []),
    ].find((option) => option.textContent?.includes(`Viridis`))
    expect(viridis_option?.querySelector(`.colorbar`)).not.toBeNull()
    viridis_option?.click()
    await tick()
    expect(new URL(location.href).searchParams.has(`color_scale`)).toBe(false)
  })

  it.each([`bogus`, `constructor`, `__proto__`])(
    `defaults an unknown filter: %s`,
    async (filter) => {
      await mount_with_url(DataTmiPage, `http://localhost/data/tmi?filter=${filter}`)

      const default_filter = radio_values().find((value) => value.startsWith(`arity=`))
      expect(checked_radio()?.value).toBe(default_filter)
    },
  )
})
