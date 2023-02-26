<script lang="ts">
  import { repository } from '$site/package.json'
  import Icon from '@iconify/svelte'
</script>

<footer>
  Questions/feedback?
  <a href="{repository}/issues"><Icon icon="octicon:mark-github" inline /></a>
  <a href="mailto:janosh@lbl.gov?subject=Matbench Discovery">
    <Icon icon="mdi:email" inline />
  </a>
</footer>

<style>
  footer {
    display: flex;
    gap: 1ex;
    place-content: center;
    place-items: center;
    padding: 3vh 3vw;
    background: #00061a;
  }
</style>
