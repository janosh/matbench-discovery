import { dirname } from 'path'
import type { PageServerLoad } from './$types'

export const load: PageServerLoad = async () => {
  const model_metas = Object.entries(
    import.meta.glob(`$root/models/**/metadata.yml`, {
      eager: true,
    })
  ).map(([key, module]) => [dirname(key), module.default])

  return { model_metas }
}
