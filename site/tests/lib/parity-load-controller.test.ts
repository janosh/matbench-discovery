import { ParityLoadController } from '$lib/parity/load-controller.svelte'
import { expect, it } from 'vitest'

it(`cancels an older request when cached data becomes ready`, async () => {
  const controller = new ParityLoadController()
  let release_load = (): void => {}
  const pending_load = controller.run(async (is_current) => {
    await new Promise<void>((resolve) => {
      release_load = resolve
    })
    expect(is_current()).toBe(false)
  })

  expect(controller.status).toBe(`loading`)
  controller.set_ready()
  release_load()
  await pending_load

  expect(controller.status).toBe(`ready`)
  expect(controller.error_message).toBe(``)
})

it(`exposes errors from the current request`, async () => {
  const controller = new ParityLoadController()
  await controller.run(async () => {
    throw new Error(`asset unavailable`)
  })

  expect(controller.status).toBe(`error`)
  expect(controller.error_message).toBe(`asset unavailable`)
})
