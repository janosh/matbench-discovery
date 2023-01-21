<script lang="ts">
  export let src: string
  export let alt: string
</script>

<figure>
  <img src={src.replace(`.`, ``)} {alt} />
  <!-- <figcaption>{alt}</figcaption> -->
</figure>

<style>
  figure {
    display: block;
  }
  figcaption {
    text-align: center;
    margin: 1ex 1em;
  }
</style>
