import schema from '$root/tests/model-metadata-schema.yml'
import fs from 'fs'
import { compile } from 'json-schema-to-typescript'

// keep model-metadata.d.ts in sync with SOT model-metadata-schema.yml
// i.e. use json-schema-to-typescript to auto-convert YAML schema to TypeScript interface
const model_md_type = await compile(schema, `ModelMetadata`, {
  style: { singleQuote: true, semi: false, printWidth: 100 },
})
const dts_out_file = `src/lib/model-metadata.d.ts`
fs.writeFileSync(dts_out_file, model_md_type)
