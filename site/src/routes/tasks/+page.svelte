<script lang="ts">
  import MODELING_TASKS from '$pkg/modeling-tasks.yml'
  import { type IconName, SubpageGrid } from 'matterviz'

  const task_icons: Record<string, IconName> = {
    discovery: `Search`,
    geo_opt: `RulerSquareCompass`,
    phonons: `Phonons`,
    diatomics: `Molecule`,
    md: `Thermometer`,
  }
  const subpages = Object.entries(MODELING_TASKS)
    .filter(([key]) => key in task_icons) // CPS has no dedicated page
    .map(([key, task]) => ({
      href: `/tasks/${key.replaceAll(`_`, `-`)}`,
      title: task.label,
      description: task.description,
      icon: task_icons[key],
    }))
</script>

<SubpageGrid
  title="Benchmark Tasks"
  subtitle="Each task probes a different aspect of machine-learning interatomic potentials, from ground-state stability prediction to finite-temperature dynamics."
  {subpages}
/>
