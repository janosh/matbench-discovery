import { writeFileSync } from 'fs'
import { dirname } from 'path'
import { fileURLToPath } from 'url'

const dir_name = dirname(fileURLToPath(import.meta.url))

for (const [path, file] of Object.entries(
  import.meta.glob(`./*.md`, { as: `raw`, eager: true })
)) {
  writeFileSync(
    `${dir_name}/${path}`,
    file
      .replaceAll(`<b>`, ``)
      .replaceAll(`</b>`, ``)
      .replaceAll(
        `src="https://img.shields.io/badge/-source-cccccc?style=flat-square"`,
        `src="https://img.shields.io/badge/source-blue?style=flat"`
      )
  )
}
