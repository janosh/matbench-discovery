import type { ObjectOption } from 'svelte-widgets'
import { untrack } from 'svelte'

type ModelSelectionConfig<T extends ObjectOption> = {
  options: T[] // canonical order, option.value = page's model identifier
  defaults: string[] // values selected when the URL has no `models` param
}

// Multi-select model choice encoded in a `models` URL query param as a comma list of
// tokens: parses the list (dropping unknown tokens), serializes selections in canonical
// option order, and reports the default serialization so bind_url_params omits the
// param when the selection matches the defaults. Config is a thunk so pages can pass
// $derived options/defaults; read() re-syncs on every navigation.
export class UrlModelSelection<T extends ObjectOption = ObjectOption> {
  selected: T[] = $state([])

  constructor(private readonly config: () => ModelSelectionConfig<T>) {
    this.selected = untrack(() => this.options_for(config().defaults))
  }

  get values(): string[] {
    return this.selected.map((option) => String(option.value))
  }

  read = (params: URLSearchParams): void => {
    const model_param = params.get(`models`)
    this.selected = this.options_for(
      model_param === null ? this.config().defaults : model_param.split(`,`),
    )
  }

  get url_entry(): [key: string, value: string, default_value: string] {
    const { defaults } = this.config()
    return [`models`, this.param_value(this.values), this.param_value(defaults)]
  }

  private options_for(values: string[]): T[] {
    return this.config().options.filter((opt) => values.includes(String(opt.value)))
  }

  private param_value(values: string[]): string {
    return this.options_for(values)
      .map((opt) => String(opt.value))
      .join(`,`)
  }
}
