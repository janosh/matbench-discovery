<script lang="ts">
  import Select from 'svelte-multiselect'
  import type { MultiSelectProps, ObjectOption, Option } from 'svelte-multiselect'

  let {
    options,
    selected = $bindable([]),
    placeholder = `Select models to plot`,
    style = `width: fit-content; max-width: 100%; min-width: min(19rem, 100%); border: 1px solid var(--border)`,
    ...rest
  }: {
    options: ObjectOption[]
    selected?: ObjectOption[]
    placeholder?: string
    style?: string
    // forward all other Select props (Option-generic so handlers stay assignable)
  } & Omit<
    MultiSelectProps<Option>,
    `options` | `selected` | `placeholder` | `style`
  > = $props()
</script>

<Select class="model-select" {options} bind:selected {placeholder} {style} {...rest} />
