<script lang="ts">
  import { onMount } from 'svelte'
  import { sortable } from 'svelte-zoo/actions'

  let table: HTMLTableElement

  // highlight all occurrences of _71_
  onMount(() => {
    const rows = table.getElementsByTagName(`tr`)
    for (const row of rows) {
      const cells = row.getElementsByTagName(`td`)
      for (const cell of cells) {
        // replace _71_ with _<em>71</em>_ to highlight
        cell.innerHTML = cell.innerText.replace(`_71_`, `_<em>71</em>_`)
      }
    }
  })
</script>

<table use:sortable {...$$props} id="T_68e58" bind:this={table}>
  <thead>
    <tr>
      <th class="blank level0">&nbsp;</th>
      <th id="T_68e58_level0_col0" class="col_heading level0 col0"
        >Isopointal Prototypes in Shaded Rectangle</th
      >
      <th id="T_68e58_level0_col1" class="col_heading level0 col1">count</th>
      <th id="T_68e58_level0_col2" class="col_heading level0 col2">MP occurrences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T_68e58_level0_row0" class="row_heading level0 row0">0</th>
      <td id="T_68e58_row0_col0" class="data row0 col0">AB2C2_oI10_71_a_f_h</td>
      <td id="T_68e58_row0_col1" class="data row0 col1">456</td>
      <td id="T_68e58_row0_col2" class="data row0 col2">75</td>
    </tr>
    <tr>
      <th id="T_68e58_level0_row1" class="row_heading level0 row1">1</th>
      <td id="T_68e58_row1_col0" class="data row1 col0">ABC2_oI8_71_a_b_f</td>
      <td id="T_68e58_row1_col1" class="data row1 col1">194</td>
      <td id="T_68e58_row1_col2" class="data row1 col2">0</td>
    </tr>
    <tr>
      <th id="T_68e58_level0_row2" class="row_heading level0 row2">2</th>
      <td id="T_68e58_row2_col0" class="data row2 col0">AB2_oI6_71_a_e</td>
      <td id="T_68e58_row2_col1" class="data row2 col1">68</td>
      <td id="T_68e58_row2_col2" class="data row2 col2">16</td>
    </tr>
    <tr>
      <th id="T_68e58_level0_row3" class="row_heading level0 row3">3</th>
      <td id="T_68e58_row3_col0" class="data row3 col0">ABC3_cP5_221_a_b_c</td>
      <td id="T_68e58_row3_col1" class="data row3 col1">30</td>
      <td id="T_68e58_row3_col2" class="data row3 col2">1429</td>
    </tr>
    <tr>
      <th id="T_68e58_level0_row4" class="row_heading level0 row4">4</th>
      <td id="T_68e58_row4_col0" class="data row4 col0">ABC2_cF16_225_a_b_c</td>
      <td id="T_68e58_row4_col1" class="data row4 col1">24</td>
      <td id="T_68e58_row4_col2" class="data row4 col2">4395</td>
    </tr>
    <tr>
      <th id="T_68e58_level0_row5" class="row_heading level0 row5">5</th>
      <td id="T_68e58_row5_col0" class="data row5 col0">AB3_oI8_71_a_bf</td>
      <td id="T_68e58_row5_col1" class="data row5 col1">19</td>
      <td id="T_68e58_row5_col2" class="data row5 col2">2</td>
    </tr>
    <tr>
      <th id="T_68e58_level0_row6" class="row_heading level0 row6">6</th>
      <td id="T_68e58_row6_col0" class="data row6 col0">AB2_mC12_15_e_2f</td>
      <td id="T_68e58_row6_col1" class="data row6 col1">18</td>
      <td id="T_68e58_row6_col2" class="data row6 col2">0</td>
    </tr>
    <tr>
      <th id="T_68e58_level0_row7" class="row_heading level0 row7">7</th>
      <td id="T_68e58_row7_col0" class="data row7 col0">ABC2_tI8_139_a_b_e</td>
      <td id="T_68e58_row7_col1" class="data row7 col1">17</td>
      <td id="T_68e58_row7_col2" class="data row7 col2">5</td>
    </tr>
    <tr>
      <th id="T_68e58_level0_row8" class="row_heading level0 row8">8</th>
      <td id="T_68e58_row8_col0" class="data row8 col0">AB2C5_oI16_71_c_e_am</td>
      <td id="T_68e58_row8_col1" class="data row8 col1">16</td>
      <td id="T_68e58_row8_col2" class="data row8 col2">3</td>
    </tr>
    <tr>
      <th id="T_68e58_level0_row9" class="row_heading level0 row9">9</th>
      <td id="T_68e58_row9_col0" class="data row9 col0">A2B3_oI10_71_e_cg</td>
      <td id="T_68e58_row9_col1" class="data row9 col1">13</td>
      <td id="T_68e58_row9_col2" class="data row9 col2">3</td>
    </tr>
  </tbody>
</table>

<style type="text/css">
  #T_68e58_row0_col1,
  #T_68e58_row4_col2 {
    background-color: #fde725;
    color: #000000;
  }
  #T_68e58_row0_col2 {
    background-color: #46075a;
    color: #f1f1f1;
  }
  #T_68e58_row1_col1 {
    background-color: #297a8e;
    color: #f1f1f1;
  }
  #T_68e58_row1_col2,
  #T_68e58_row2_col2,
  #T_68e58_row5_col2,
  #T_68e58_row6_col2,
  #T_68e58_row7_col2,
  #T_68e58_row8_col2,
  #T_68e58_row9_col1,
  #T_68e58_row9_col2 {
    background-color: #440154;
    color: #f1f1f1;
  }
  #T_68e58_row2_col1 {
    background-color: #472c7a;
    color: #f1f1f1;
  }
  #T_68e58_row3_col1 {
    background-color: #470e61;
    color: #f1f1f1;
  }
  #T_68e58_row3_col2 {
    background-color: #31668e;
    color: #f1f1f1;
  }
  #T_68e58_row4_col1 {
    background-color: #460a5d;
    color: #f1f1f1;
  }
  #T_68e58_row5_col1 {
    background-color: #450559;
    color: #f1f1f1;
  }
  #T_68e58_row6_col1,
  #T_68e58_row7_col1 {
    background-color: #450457;
    color: #f1f1f1;
  }
  #T_68e58_row8_col1 {
    background-color: #440256;
    color: #f1f1f1;
  }
  table :global(em) {
    font-weight: bolder;
    color: orange;
    font-style: normal;
  }
</style>
