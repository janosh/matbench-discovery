// Parse YAML environment.dependencies specs into displayable package links.

export type ParsedDependency = { name: string; detail: string; href: string }

// PEP 508 name, optionally with extras like pkg[extra]
const PKG = `([A-Za-z0-9][A-Za-z0-9._-]*(?:\\[[^\\]]+\\])?)`
const AT_SPEC = new RegExp(`^${PKG}\\s*@\\s*(.+)$`)
const PIN_SPEC = new RegExp(`^${PKG}\\s*(==|>=|<=|!=|~=|<|>)\\s*(.+)$`)

const pypi_href = (name: string, version = ``): string =>
  `https://pypi.org/project/${name.replace(/\[[^\]]*\]$/, ``)}/${version}`

export function parse_dependency_spec(dep: string): ParsedDependency {
  const trimmed = dep.trim()
  const at_match = AT_SPEC.exec(trimmed)
  if (at_match) {
    const name = at_match[1]
    const locator = at_match[2].trim()
    const url = locator.replace(/^git\+/, ``) // "git+https://..." -> "https://..."
    const href = url.startsWith(`http`) ? url : pypi_href(name)
    return { name, detail: locator, href }
  }
  const pin_match = PIN_SPEC.exec(trimmed)
  if (pin_match) {
    const name = pin_match[1]
    const operator = pin_match[2]
    const version = pin_match[3].trim()
    const href = pypi_href(name, operator === `==` ? version : undefined)
    return { name, detail: `${operator}${version}`, href }
  }
  return { name: trimmed, detail: ``, href: pypi_href(trimmed) }
}
