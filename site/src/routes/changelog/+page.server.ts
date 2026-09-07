import changelog from '$root/changelog.md?raw'
import { assert_ok, create_markdown, render_markdown } from 'svelte-widgets/markdown'

const decrease_heading_level = (str: string) => str.replaceAll(`###`, `#`)

export const load = async () => {
  const engine = create_markdown()
  const document = assert_ok(
    await engine.parse(decrease_heading_level(changelog), {
      filename: `changelog.md`,
      dialect: `markdown`,
    }),
  )
  return { changelog: assert_ok(await render_markdown(document)) }
}
