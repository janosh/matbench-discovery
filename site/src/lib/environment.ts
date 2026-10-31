// Parse YAML environment.dependencies specs into displayable package links.

// PEP 508 name, optionally with extras like pkg[extra]
const PKG = `([A-Za-z0-9][A-Za-z0-9._-]*(?:\\[[^\\]]+\\])?)`
const AT_SPEC = new RegExp(`^${PKG}\\s*@\\s*(.+)$`)
const PIN_SPEC = new RegExp(`^${PKG}\\s*(==|>=|<=|!=|~=|<|>)\\s*(.+)$`)

// local version labels (torch's +cu126) have no PyPI release page
const pypi_href = (name: string, version = ``): string =>
  `https://pypi.org/project/${name.replace(/\[[^\]]*\]$/, ``)}/${version.replace(/\+.*$/, ``)}`

const vcs_browser_href = (locator: string): string => {
  const url = locator.replace(/^git\+/, ``)
  const github_match =
    /^https:\/\/github\.com\/(?<owner>[^/]+)\/(?<repo>[^/@]+)@(?<revision>[^#]+)(?<fragment>#.*)?$/.exec(
      url,
    )
  if (!github_match?.groups) return url

  const { owner, revision, fragment = `` } = github_match.groups
  const repo = github_match.groups.repo.replace(/\.git$/, ``)
  const revision_kind = /^[0-9a-f]{7,40}$/i.test(revision) ? `commit` : `tree`
  return `https://github.com/${owner}/${repo}/${revision_kind}/${revision}${fragment}`
}

export function parse_dependency_spec(dep: string) {
  const trimmed = dep.trim()
  const at_match = AT_SPEC.exec(trimmed)
  if (at_match) {
    const name = at_match[1]
    const locator = at_match[2].trim()
    const url = vcs_browser_href(locator)
    const href = url.startsWith(`http`) ? url : pypi_href(name)
    // a VCS install's name links its repo (may not be on PyPI), minus any user@ credentials
    const repo_match =
      /^git\+(?<scheme>https?:\/\/)(?:[^@/]+@)?(?<path>[^@]+?)(?:\.git)?@/.exec(locator)
    const repo =
      repo_match?.groups && `${repo_match.groups.scheme}${repo_match.groups.path}`
    return { name, name_href: repo ?? pypi_href(name), detail: locator, href }
  }
  const pin_match = PIN_SPEC.exec(trimmed)
  if (pin_match) {
    const name = pin_match[1]
    const operator = pin_match[2]
    const version = pin_match[3].trim()
    const href = pypi_href(name, operator === `==` ? version : undefined)
    return { name, name_href: pypi_href(name), detail: `${operator}${version}`, href }
  }
  const href = pypi_href(trimmed)
  return { name: trimmed, name_href: href, detail: ``, href }
}
