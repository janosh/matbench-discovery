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
