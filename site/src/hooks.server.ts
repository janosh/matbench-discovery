import { preprint } from '$site/package.json'
import type { Handle } from '@sveltejs/kit'
import { redirect } from '@sveltejs/kit'

export const handle: Handle = async ({ event, resolve }) => {
  // Redirect /preprint to the arXiv paper
  if (event.url.pathname === `/preprint`) {
    throw redirect(307, preprint)
  }

  const response = await resolve(event)
  return response
}
