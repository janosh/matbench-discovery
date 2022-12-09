import { redirect } from '@sveltejs/kit'
import type { LayoutServerLoad } from './$types'

export const load: LayoutServerLoad = ({ url }) => {
  if (url.pathname === '/api') {
    throw redirect(307, '/api/data')
  }
}
