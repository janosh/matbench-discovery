<script lang="ts">
  import MODELING_TASKS from '$pkg/modeling-tasks.yml'
  import { SubpageGrid } from 'matterviz'
  import {
    Molecule,
    Phonons,
    RulerSquareCompass,
    Search,
    Thermometer,
  } from 'svelte-widgets/icons'

  // CPS has no dedicated page
  const task_page_icons = [
    [`discovery`, Search],
    [`geo_opt`, RulerSquareCompass],
    [`phonons`, Phonons],
    [`diatomics`, Molecule],
    [`md`, Thermometer],
  ] as const
  const subpages = task_page_icons.map(([key, icon]) => ({
    href: `/tasks/${key.replaceAll(`_`, `-`)}`,
    title: MODELING_TASKS[key].label,
    description: MODELING_TASKS[key].description,
    icon,
  }))
</script>

<SubpageGrid
  title="Benchmark Tasks"
  subtitle="Each task probes a different aspect of ML force fields, from ground-state stability prediction to finite-temperature dynamics."
  {subpages}
/>
