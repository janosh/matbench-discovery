<script lang="ts">
  import { Icon } from '$lib'
  import pkg from '$site/package.json'
</script>

<footer>
  <nav>
    <a href="{pkg.repository}/issues">
      <Icon icon="Issues" />
      Issues
    </a>
    <a href="mailto:janosh.riebesell@gmail.com?subject=Matbench Discovery">
      <Icon icon="Contact" />
      Contact
    </a>
    <a href="/changelog">
      <Icon icon="Changelog" />
      Changelog
    </a>
    <a href="/rss.xml" title="Be notified of new model submissions">
      <Icon icon="RSS" />
      RSS
    </a>
  </nav>
  <img src="/favicon.svg" alt="Logo" width="30px" style="vertical-align: middle" />
  &ensp;{pkg.title} &ensp; | &ensp; ©
  <a href={pkg[`author-url`]}>{pkg.author.split(`<`)[0]}</a>
  (<a href="{pkg.repository}/blob/-/license">2022</a>)
</footer>

<style>
  footer {
    padding: 3vh 3vw;
    background: #00061a;
    text-align: center;
  }
  footer nav {
    display: flex;
    gap: 2em;
    justify-content: center;
    margin: 2em 0;
  }
  footer > nav svg {
    margin-right: 3pt;
  }
</style>
