<script lang="ts">
  import Select from 'svelte-multiselect'
  import type { ObjectOption } from 'svelte-multiselect'

  let {
    options,
    selected = $bindable([]),
    placeholder = `Select models to plot`,
  }: {
    options: ObjectOption[]
    selected?: ObjectOption[]
    placeholder?: string
  } = $props()
</script>

<Select
  class="model-select"
  {options}
  bind:selected
  {placeholder}
  style="width: fit-content; max-width: 100%; min-width: min(19rem, 100%); border: 1px solid var(--border)"
/>
