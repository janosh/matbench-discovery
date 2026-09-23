// Measure the unwrapped toolbar so optional links cannot push primary controls
// onto another line. Ignore the filter chips, which have their own row.
export function fit_toolbar_links(controls: HTMLElement) {
  const toolbar = controls.parentElement
  if (!toolbar?.matches(`.control-buttons`)) return undefined
  const links = [...toolbar.querySelectorAll<HTMLElement>(`[data-toolbar-optional]`)]
  if (!links.length) return undefined

  const update = () => {
    for (const link of links) link.hidden = false
    const available_width = toolbar.clientWidth
    const chips = toolbar.querySelector<HTMLElement>(`:scope > .active-filters`)
    const styles = [toolbar, controls, ...(chips ? [chips] : [])].map(
      (element) => [element, element.style.cssText] as const,
    )
    // Measure out of flow so the toolbar cannot widen its table while measuring.
    Object.assign(toolbar.style, {
      position: `absolute`,
      inlineSize: `max-content`,
      contain: `none`,
      flexWrap: `nowrap`,
    })
    Object.assign(controls.style, { flex: `0 0 auto`, flexWrap: `nowrap` })
    if (chips) chips.style.display = `none`
    const required_width = toolbar.getBoundingClientRect().width
    for (const [element, style] of styles) element.style.cssText = style
    for (const link of links) link.hidden = required_width > available_width
  }
  const resize = new ResizeObserver(update)
  resize.observe(toolbar)
  const mutations = new MutationObserver(update)
  mutations.observe(controls, { childList: true, characterData: true, subtree: true })
  update()
  return () => {
    resize.disconnect()
    mutations.disconnect()
  }
}
