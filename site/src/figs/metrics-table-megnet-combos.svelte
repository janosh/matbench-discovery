<script lang="ts">
  import { onMount } from 'svelte'
  import { sortable } from 'svelte-zoo/actions'

  export let hide: string[] = []

  onMount(() => {
    const rows = document.querySelectorAll(`tr`)
    for (const row of rows) {
      const header = row.children[0].textContent
      if (header && hide.includes(header)) {
        row.style.display = `none`
      }
    }
  })
</script>

<table use:sortable {...$$props} id="T_">
  <thead>
    <tr>
      <th class="index_name level0">model</th>
      <th id="T__level0_col0" class="col_heading level0 col0">F1</th>
      <th id="T__level0_col1" class="col_heading level0 col1">R²</th>
      <th id="T__level0_col2" class="col_heading level0 col2">DAF</th>
      <th id="T__level0_col3" class="col_heading level0 col3">Precision</th>
      <th id="T__level0_col5" class="col_heading level0 col5">TPR</th>
      <th id="T__level0_col7" class="col_heading level0 col7">TNR</th>
      <th id="T__level0_col9" class="col_heading level0 col9">Accuracy</th>
      <th id="T__level0_col10" class="col_heading level0 col10">MAE</th>
      <th id="T__level0_col11" class="col_heading level0 col11">RMSE</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T__level0_row0" class="row_heading level0 row0">CHGNet</th>
      <td id="T__row0_col0" class="data row0 col0">0.59</td>
      <td id="T__row0_col1" class="data row0 col1">0.61</td>
      <td id="T__row0_col2" class="data row0 col2">3.06</td>
      <td id="T__row0_col3" class="data row0 col3">0.52</td>
      <td id="T__row0_col5" class="data row0 col5">0.67</td>
      <td id="T__row0_col7" class="data row0 col7">0.87</td>
      <td id="T__row0_col9" class="data row0 col9">0.84</td>
      <td id="T__row0_col10" class="data row0 col10">0.07</td>
      <td id="T__row0_col11" class="data row0 col11">0.11</td>
    </tr>
    <tr>
      <th id="T__level0_row1" class="row_heading level0 row1">M3GNet</th>
      <td id="T__row1_col0" class="data row1 col0">0.58</td>
      <td id="T__row1_col1" class="data row1 col1">0.59</td>
      <td id="T__row1_col2" class="data row1 col2">2.66</td>
      <td id="T__row1_col3" class="data row1 col3">0.45</td>
      <td id="T__row1_col5" class="data row1 col5">0.79</td>
      <td id="T__row1_col7" class="data row1 col7">0.80</td>
      <td id="T__row1_col9" class="data row1 col9">0.80</td>
      <td id="T__row1_col10" class="data row1 col10">0.07</td>
      <td id="T__row1_col11" class="data row1 col11">0.12</td>
    </tr>
    <tr>
      <th id="T__level0_row2" class="row_heading level0 row2">M3GNet + MEGNet</th>
      <td id="T__row2_col0" class="data row2 col0">0.55</td>
      <td id="T__row2_col1" class="data row2 col1">0.35</td>
      <td id="T__row2_col2" class="data row2 col2">2.54</td>
      <td id="T__row2_col3" class="data row2 col3">0.43</td>
      <td id="T__row2_col5" class="data row2 col5">0.77</td>
      <td id="T__row2_col7" class="data row2 col7">0.79</td>
      <td id="T__row2_col9" class="data row2 col9">0.79</td>
      <td id="T__row2_col10" class="data row2 col10">0.09</td>
      <td id="T__row2_col11" class="data row2 col11">0.14</td>
    </tr>
    <tr>
      <th id="T__level0_row3" class="row_heading level0 row3">CHGNet + MEGNet</th>
      <td id="T__row3_col0" class="data row3 col0">0.55</td>
      <td id="T__row3_col1" class="data row3 col1">0.34</td>
      <td id="T__row3_col2" class="data row3 col2">2.56</td>
      <td id="T__row3_col3" class="data row3 col3">0.43</td>
      <td id="T__row3_col5" class="data row3 col5">0.75</td>
      <td id="T__row3_col7" class="data row3 col7">0.79</td>
      <td id="T__row3_col9" class="data row3 col9">0.79</td>
      <td id="T__row3_col10" class="data row3 col10">0.09</td>
      <td id="T__row3_col11" class="data row3 col11">0.15</td>
    </tr>
    <tr>
      <th id="T__level0_row4" class="row_heading level0 row4">MEGNet</th>
      <td id="T__row4_col0" class="data row4 col0">0.52</td>
      <td id="T__row4_col1" class="data row4 col1">-0.27</td>
      <td id="T__row4_col2" class="data row4 col2">2.70</td>
      <td id="T__row4_col3" class="data row4 col3">0.46</td>
      <td id="T__row4_col5" class="data row4 col5">0.59</td>
      <td id="T__row4_col7" class="data row4 col7">0.86</td>
      <td id="T__row4_col9" class="data row4 col9">0.81</td>
      <td id="T__row4_col10" class="data row4 col10">0.13</td>
      <td id="T__row4_col11" class="data row4 col11">0.20</td>
    </tr>
    <tr>
      <th id="T__level0_row5" class="row_heading level0 row5">CGCNN</th>
      <td id="T__row5_col0" class="data row5 col0">0.52</td>
      <td id="T__row5_col1" class="data row5 col1">-0.61</td>
      <td id="T__row5_col2" class="data row5 col2">2.62</td>
      <td id="T__row5_col3" class="data row5 col3">0.45</td>
      <td id="T__row5_col5" class="data row5 col5">0.60</td>
      <td id="T__row5_col7" class="data row5 col7">0.85</td>
      <td id="T__row5_col9" class="data row5 col9">0.81</td>
      <td id="T__row5_col10" class="data row5 col10">0.14</td>
      <td id="T__row5_col11" class="data row5 col11">0.23</td>
    </tr>
    <tr>
      <th id="T__level0_row6" class="row_heading level0 row6">CGCNN+P</th>
      <td id="T__row6_col0" class="data row6 col0">0.51</td>
      <td id="T__row6_col1" class="data row6 col1">0.02</td>
      <td id="T__row6_col2" class="data row6 col2">2.38</td>
      <td id="T__row6_col3" class="data row6 col3">0.41</td>
      <td id="T__row6_col5" class="data row6 col5">0.69</td>
      <td id="T__row6_col7" class="data row6 col7">0.79</td>
      <td id="T__row6_col9" class="data row6 col9">0.78</td>
      <td id="T__row6_col10" class="data row6 col10">0.11</td>
      <td id="T__row6_col11" class="data row6 col11">0.18</td>
    </tr>
    <tr>
      <th id="T__level0_row7" class="row_heading level0 row7">Wrenformer</th>
      <td id="T__row7_col0" class="data row7 col0">0.48</td>
      <td id="T__row7_col1" class="data row7 col1">-0.04</td>
      <td id="T__row7_col2" class="data row7 col2">2.13</td>
      <td id="T__row7_col3" class="data row7 col3">0.36</td>
      <td id="T__row7_col5" class="data row7 col5">0.71</td>
      <td id="T__row7_col7" class="data row7 col7">0.74</td>
      <td id="T__row7_col9" class="data row7 col9">0.74</td>
      <td id="T__row7_col10" class="data row7 col10">0.10</td>
      <td id="T__row7_col11" class="data row7 col11">0.18</td>
    </tr>
    <tr>
      <th id="T__level0_row8" class="row_heading level0 row8">BOWSR + MEGNet</th>
      <td id="T__row8_col0" class="data row8 col0">0.44</td>
      <td id="T__row8_col1" class="data row8 col1">0.15</td>
      <td id="T__row8_col2" class="data row8 col2">1.90</td>
      <td id="T__row8_col3" class="data row8 col3">0.32</td>
      <td id="T__row8_col5" class="data row8 col5">0.74</td>
      <td id="T__row8_col7" class="data row8 col7">0.67</td>
      <td id="T__row8_col9" class="data row8 col9">0.68</td>
      <td id="T__row8_col10" class="data row8 col10">0.11</td>
      <td id="T__row8_col11" class="data row8 col11">0.16</td>
    </tr>
    <tr>
      <th id="T__level0_row9" class="row_heading level0 row9">Voronoi RF</th>
      <td id="T__row9_col0" class="data row9 col0">0.34</td>
      <td id="T__row9_col1" class="data row9 col1">-0.32</td>
      <td id="T__row9_col2" class="data row9 col2">1.51</td>
      <td id="T__row9_col3" class="data row9 col3">0.26</td>
      <td id="T__row9_col5" class="data row9 col5">0.52</td>
      <td id="T__row9_col7" class="data row9 col7">0.69</td>
      <td id="T__row9_col9" class="data row9 col9">0.66</td>
      <td id="T__row9_col10" class="data row9 col10">0.14</td>
      <td id="T__row9_col11" class="data row9 col11">0.21</td>
    </tr>
  </tbody>
</table>

<style type="text/css">
  #T_ {
    font-family: sans-serif;
    border-collapse: collapse;
  }
  #T_ td {
    border: none;
    padding: 4px 6px;
    white-space: nowrap;
  }
  #T_ th {
    border: none;
    padding: 4px 6px;
    white-space: nowrap;
  }
  #T_ th {
    border: 1px solid;
    border-width: 1px 0;
    text-align: left;
  }
  #T__row0_col0,
  #T__row0_col1,
  #T__row0_col2,
  #T__row0_col3,
  #T__row0_col7,
  #T__row0_col9,
  #T__row0_col10,
  #T__row0_col11,
  #T__row1_col5 {
    background-color: #fde725;
    color: #000000;
  }
  #T__row0_col5,
  #T__row6_col2,
  #T__row7_col0 {
    background-color: #1fa188;
    color: #f1f1f1;
  }
  #T__row1_col0 {
    background-color: #dfe318;
    color: #000000;
  }
  #T__row1_col1 {
    background-color: #f8e621;
    color: #000000;
  }
  #T__row1_col2 {
    background-color: #5cc863;
    color: #000000;
  }
  #T__row1_col3 {
    background-color: #56c667;
    color: #000000;
  }
  #T__row1_col7 {
    background-color: #32b67a;
    color: #f1f1f1;
  }
  #T__row1_col9 {
    background-color: #73d056;
    color: #000000;
  }
  #T__row1_col10 {
    background-color: #f6e620;
    color: #000000;
  }
  #T__row1_col11 {
    background-color: #efe51c;
    color: #000000;
  }
  #T__row2_col0,
  #T__row4_col9 {
    background-color: #a0da39;
    color: #000000;
  }
  #T__row2_col1 {
    background-color: #75d054;
    color: #000000;
  }
  #T__row2_col2,
  #T__row3_col3 {
    background-color: #35b779;
    color: #f1f1f1;
  }
  #T__row2_col3 {
    background-color: #31b57b;
    color: #f1f1f1;
  }
  #T__row2_col5 {
    background-color: #c8e020;
    color: #000000;
  }
  #T__row2_col7 {
    background-color: #22a785;
    color: #f1f1f1;
  }
  #T__row2_col9,
  #T__row3_col11,
  #T__row7_col5 {
    background-color: #48c16e;
    color: #f1f1f1;
  }
  #T__row2_col10 {
    background-color: #54c568;
    color: #000000;
  }
  #T__row2_col11,
  #T__row3_col10 {
    background-color: #4cc26c;
    color: #000000;
  }
  #T__row3_col0 {
    background-color: #9bd93c;
    color: #000000;
  }
  #T__row3_col1 {
    background-color: #6ece58;
    color: #000000;
  }
  #T__row3_col2 {
    background-color: #3aba76;
    color: #f1f1f1;
  }
  #T__row3_col5 {
    background-color: #a5db36;
    color: #000000;
  }
  #T__row3_col7,
  #T__row6_col7 {
    background-color: #25ac82;
    color: #f1f1f1;
  }
  #T__row3_col9,
  #T__row4_col0 {
    background-color: #4ac16d;
    color: #000000;
  }
  #T__row4_col1 {
    background-color: #375b8d;
    color: #f1f1f1;
  }
  #T__row4_col2 {
    background-color: #69cd5b;
    color: #000000;
  }
  #T__row4_col3 {
    background-color: #6ccd5a;
    color: #000000;
  }
  #T__row4_col5 {
    background-color: #3b518b;
    color: #f1f1f1;
  }
  #T__row4_col7 {
    background-color: #d8e219;
    color: #000000;
  }
  #T__row4_col10 {
    background-color: #414287;
    color: #f1f1f1;
  }
  #T__row4_col11 {
    background-color: #3e4989;
    color: #f1f1f1;
  }
  #T__row5_col0 {
    background-color: #44bf70;
    color: #f1f1f1;
  }
  #T__row5_col1,
  #T__row5_col11,
  #T__row8_col7,
  #T__row9_col0,
  #T__row9_col2,
  #T__row9_col3,
  #T__row9_col5,
  #T__row9_col9,
  #T__row9_col10 {
    background-color: #440154;
    color: #f1f1f1;
  }
  #T__row5_col2 {
    background-color: #4ec36b;
    color: #000000;
  }
  #T__row5_col3 {
    background-color: #50c46a;
    color: #000000;
  }
  #T__row5_col5 {
    background-color: #34618d;
    color: #f1f1f1;
  }
  #T__row5_col7 {
    background-color: #b0dd2f;
    color: #000000;
  }
  #T__row5_col9 {
    background-color: #86d549;
    color: #000000;
  }
  #T__row5_col10 {
    background-color: #482475;
    color: #f1f1f1;
  }
  #T__row6_col0 {
    background-color: #3fbc73;
    color: #f1f1f1;
  }
  #T__row6_col1,
  #T__row7_col10 {
    background-color: #1f958b;
    color: #f1f1f1;
  }
  #T__row6_col3 {
    background-color: #1fa187;
    color: #f1f1f1;
  }
  #T__row6_col5 {
    background-color: #27ad81;
    color: #f1f1f1;
  }
  #T__row6_col9 {
    background-color: #2fb47c;
    color: #f1f1f1;
  }
  #T__row6_col10 {
    background-color: #24878e;
    color: #f1f1f1;
  }
  #T__row6_col11 {
    background-color: #277f8e;
    color: #f1f1f1;
  }
  #T__row7_col1 {
    background-color: #23898e;
    color: #f1f1f1;
  }
  #T__row7_col2,
  #T__row8_col0 {
    background-color: #2a788e;
    color: #f1f1f1;
  }
  #T__row7_col3 {
    background-color: #29798e;
    color: #f1f1f1;
  }
  #T__row7_col7 {
    background-color: #2d718e;
    color: #f1f1f1;
  }
  #T__row7_col9 {
    background-color: #26818e;
    color: #f1f1f1;
  }
  #T__row7_col11 {
    background-color: #2b748e;
    color: #f1f1f1;
  }
  #T__row8_col1 {
    background-color: #28ae80;
    color: #f1f1f1;
  }
  #T__row8_col2 {
    background-color: #3a538b;
    color: #f1f1f1;
  }
  #T__row8_col3 {
    background-color: #3f4788;
    color: #f1f1f1;
  }
  #T__row8_col5 {
    background-color: #8bd646;
    color: #000000;
  }
  #T__row8_col9 {
    background-color: #482979;
    color: #f1f1f1;
  }
  #T__row8_col10 {
    background-color: #2a768e;
    color: #f1f1f1;
  }
  #T__row8_col11 {
    background-color: #1fa088;
    color: #f1f1f1;
  }
  #T__row9_col1 {
    background-color: #3c508b;
    color: #f1f1f1;
  }
  #T__row9_col7 {
    background-color: #482878;
    color: #f1f1f1;
  }
  #T__row9_col11 {
    background-color: #424086;
    color: #f1f1f1;
  }
</style>
