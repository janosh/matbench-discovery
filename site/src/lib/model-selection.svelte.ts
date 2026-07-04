import type { ObjectOption } from 'svelte-multiselect'
import { untrack } from 'svelte'

type ModelSelectionConfig = {
  options: ObjectOption[] // canonical order, option.value = page's model identifier
  defaults: string[] // values selected when the URL has no `models` param
  // URL token -> value; return undefined to drop unknown tokens (default: identity)
  from_url?: (token: string) => string | undefined
  // value -> URL token, e.g. display name -> model key (default: identity)
  to_url?: (value: string) => string
}

// Multi-select model choice encoded in a `models` URL query param as a comma list of
// tokens: parses the list (dropping unknown tokens), serializes selections in canonical
// option order, and reports the default serialization so bind_url_params omits the
// param when the selection matches the defaults. Config is a thunk so pages can pass
// $derived options/defaults; read() re-syncs on every navigation.
export class UrlModelSelection {
  selected: ObjectOption[] = $state([])

  constructor(private config: () => ModelSelectionConfig) {
    this.selected = untrack(() => this.options_for(config().defaults))
  }

  get values(): string[] {
    return this.selected.map((option) => String(option.value))
  }

  read = (params: URLSearchParams): void => {
    const { defaults, from_url = (token: string) => token } = this.config()
    const model_param = params.get(`models`)
    const values = model_param === null ? defaults : model_param
      .split(`,`)
      .map(from_url)
      .filter((value): value is string => value !== undefined)
    this.selected = this.options_for(values)
  }

  get url_entry(): [key: string, value: string, default_value: string] {
    const { defaults } = this.config()
    return [`models`, this.param_value(this.values), this.param_value(defaults)]
  }

  private options_for(values: string[]): ObjectOption[] {
    return this.config().options.filter((opt) => values.includes(String(opt.value)))
  }

  private param_value(values: string[]): string {
    const { options, to_url = (value: string) => value } = this.config()
    return options
      .map((opt) => String(opt.value))
      .filter((value) => values.includes(value))
      .map(to_url)
      .join(`,`)
  }
}
