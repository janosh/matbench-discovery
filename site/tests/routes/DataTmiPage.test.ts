import DataTmiPage from '$routes/data/tmi/+page.svelte'
import { describe, expect, it } from 'vitest'
import { checkbox_for, mount_with_url } from '../index'

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
  })

  it(`restores filter and toggles from URL params`, async () => {
    await mount_with_url(DataTmiPage, `http://localhost/data/tmi`)
    // pick a non-default filter value from the rendered radios, then remount with it
    const target_filter = radio_values().at(-1)
    if (!target_filter) throw new Error(`no filter radios rendered`)
    document.body.innerHTML = ``

    await mount_with_url(
      DataTmiPage,
      `http://localhost/data/tmi?filter=${encodeURIComponent(
        target_filter,
      )}&normalized=1&log=1`,
    )

    expect(checked_radio()?.value).toBe(target_filter)
    expect(checkbox_for(`Normalize by data set size`).checked).toBe(true)
  })

  it(`falls back to the default filter for unknown values`, async () => {
    await mount_with_url(DataTmiPage, `http://localhost/data/tmi?filter=bogus`)

    const default_filter = radio_values().find((value) => value.startsWith(`arity=`))
    expect(checked_radio()?.value).toBe(default_filter)
  })
})
