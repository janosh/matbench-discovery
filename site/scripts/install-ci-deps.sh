set -euo pipefail

corepack enable
# Resolve main once: pnpm authorizes Git builds by source URL, not package name.
printf 'allowBuilds:\n' > pnpm-workspace.yaml
for package in matterviz svelte-widgets; do
  commit=$(git ls-remote "https://github.com/janosh/$package.git" refs/heads/main | cut -f1)
  [[ $commit =~ ^[0-9a-f]{40}$ ]]
  npm pkg set "devDependencies.$package=github:janosh/$package#$commit"
  printf '  %s@https://codeload.github.com/janosh/%s/tar.gz/%s: true\n' "$package" "$package" "$commit" >> pnpm-workspace.yaml
done
# MatterViz's Git subdependency must use the same approved widgets revision.
printf 'overrides:\n  svelte-widgets: %s\n' "$(npm pkg get devDependencies.svelte-widgets)" >> pnpm-workspace.yaml
pnpm install --config.strict-dep-builds=false --config.block-exotic-subdeps=false
