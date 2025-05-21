<script lang="ts">
  import pkg from '$site/package.json'
</script>

<footer>
  <nav>
    <a href="{pkg.repository}/issues">
      <svg><use href="#icon-issues" /></svg>
      Issues
    </a>
    <a href="mailto:janosh.riebesell@gmail.com?subject=Matbench Discovery">
      <svg><use href="#icon-contact" /></svg>
      Contact
    </a>
    <a href="/changelog">
      <svg><use href="#icon-changelog" /></svg>
      Changelog
    </a>
    <a href="/rss.xml" title="Be notified of new model submissions">
      <svg><use href="#icon-rss" /></svg>
      RSS
    </a>
  </nav>
  <img src="/favicon.svg" alt="Logo" width="30px" style="vertical-align: middle;" />
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
