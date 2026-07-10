import type { SvelteSet } from 'svelte/reactivity'

// Attachment factory: track which plot keys are near the viewport so pages with
// dozens of per-element plots only render the visible ones. Falls back to
// always-visible when IntersectionObserver is unavailable (e.g. happy-dom tests).
export const make_plot_observer =
  (visible: SvelteSet<string>) => (key: string) => (node: HTMLElement) => {
    const observer = globalThis.IntersectionObserver
      ? new IntersectionObserver(
          ([entry]) => {
            if (entry) visible[entry.isIntersecting ? `add` : `delete`](key)
          },
          { rootMargin: `300px 0px` },
        )
      : null
    if (observer) observer.observe(node)
    else visible.add(key)
    return () => {
      visible.delete(key)
      observer?.disconnect()
    }
  }
