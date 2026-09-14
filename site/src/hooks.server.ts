import { building } from '$app/environment'
import pkg from '$site/package.json'
import type { Handle } from '@sveltejs/kit'
import { redirect } from '@sveltejs/kit'

export const handle: Handle = async ({ event, resolve }) => {
  if (event.url.pathname === `/preprint`) {
    redirect(307, pkg.preprint)
  }
  if (event.url.pathname === `/tasks`) {
    redirect(307, `/benchmarks${building ? `` : event.url.search}`)
  }

  return resolve(event)
}
