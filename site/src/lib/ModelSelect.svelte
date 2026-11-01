<script lang="ts" module>
  import { pick_contrast_color } from 'svelte-widgets/attachments'

  // Selected-chip style in a model's color. Text turns white once white reaches 3:1, as
  // WCAG's pure max-ratio pick puts black on mid-tones (e.g. #2ca02c) where white reads better.
  export const model_chip_style = (color: string) => ({
    selected: `background: ${color}; color: ${pick_contrast_color({ bg_color: color, min_ratio: 3 })};`,
    option: ``,
  })
</script>

<script lang="ts" generics="T extends Option">
  import { MultiSelect } from 'svelte-widgets'
  import type { MultiSelectProps, Option } from 'svelte-widgets'

  let {
    options,
    value = $bindable([]),
    placeholder = `Select models to plot`,
    style = `width: fit-content; max-width: min(48rem, 100%); min-width: min(19rem, 100%); border: 1px solid var(--border)`,
    ...rest
  }: Extract<MultiSelectProps<T>, { mode?: `multiple` }> = $props()
</script>

<MultiSelect {options} bind:value {placeholder} {style} {...rest} />
