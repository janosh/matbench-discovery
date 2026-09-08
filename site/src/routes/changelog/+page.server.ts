import changelog from '$root/changelog.md?raw'
import { assert_ok, create_markdown } from 'svelte-widgets/markdown'

const decrease_heading_level = (str: string) => str.replaceAll(`###`, `#`)

export const load = async () => ({
  changelog: assert_ok(
    await create_markdown().render(decrease_heading_level(changelog), {
      filename: `changelog.md`,
    }),
  ),
})
