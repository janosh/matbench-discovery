<script lang="ts" generics="T extends Option">
  import Select from 'svelte-multiselect'
  import type { MultiSelectProps, Option } from 'svelte-multiselect'

  // only the bound props need the precise T; pass-through props (handlers etc.) use the
  // base Option generic so spreading them onto Select type-checks without a cast
  let {
    options,
    selected = $bindable([]),
    value = $bindable(null),
    placeholder = `Select models to plot`,
    style = `width: fit-content; max-width: min(48rem, 100%); min-width: min(19rem, 100%); border: 1px solid var(--border)`,
    ...rest
  }: Pick<MultiSelectProps<T>, `options` | `selected` | `value`> &
    Omit<MultiSelectProps<Option>, `options` | `selected` | `value`> = $props()
</script>

<Select {options} bind:selected bind:value {placeholder} {style} {...rest} />
