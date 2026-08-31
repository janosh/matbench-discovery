<script lang="ts">
  import MODELING_TASKS from '$pkg/modeling-tasks.yml'
  import type { ComponentProps } from 'svelte'
  import { SubpageGrid } from 'svelte-widgets'
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
  // SubpageGrid takes [title, href, description, icon] tuples
  const subpages: ComponentProps<typeof SubpageGrid>[`subpages`] = task_page_icons.map(
    ([key, icon]) => [
      MODELING_TASKS[key].label,
      `/tasks/${key.replaceAll(`_`, `-`)}`,
      MODELING_TASKS[key].description,
      icon,
    ],
  )
</script>

<SubpageGrid
  title="Benchmark Tasks"
  subtitle="Each task probes a different aspect of ML force fields, from ground-state stability prediction to finite-temperature dynamics."
  {subpages}
/>
