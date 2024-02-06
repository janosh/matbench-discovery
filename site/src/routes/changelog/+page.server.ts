import changelog from '$root/changelog.md?raw'
import { compile } from 'mdsvex'

const decrease_heading_level = (str: string) => str.replaceAll(`###`, `#`)

export const load = async () => ({
  changelog: await compile(decrease_heading_level(changelog)),
})
