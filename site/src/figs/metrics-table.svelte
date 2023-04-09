<script lang="ts">
  import { sortable } from 'svelte-zoo/actions'
</script>

<table use:sortable class="roomy" {...$$props} id="T_">
  <thead>
    <tr>
      <th class="blank level0">&nbsp;</th>
      <th id="T__level0_col0" class="col_heading level0 col0">F1</th>
      <th id="T__level0_col1" class="col_heading level0 col1">DAF</th>
      <th id="T__level0_col2" class="col_heading level0 col2">Precision</th>
      <th id="T__level0_col4" class="col_heading level0 col4">TPR</th>
      <th id="T__level0_col6" class="col_heading level0 col6">TNR</th>
      <th id="T__level0_col8" class="col_heading level0 col8">Accuracy</th>
      <th id="T__level0_col9" class="col_heading level0 col9">MAE</th>
      <th id="T__level0_col10" class="col_heading level0 col10">RMSE</th>
      <th id="T__level0_col11" class="col_heading level0 col11">R<sup>2</sup></th>
      <th id="T__level0_col14" class="col_heading level0 col14">model class</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T__level0_row0" class="row_heading level0 row0">CHGNet</th>
      <td id="T__row0_col0" class="data row0 col0">0.59</td>
      <td id="T__row0_col1" class="data row0 col1">3.06</td>
      <td id="T__row0_col2" class="data row0 col2">0.52</td>
      <td id="T__row0_col4" class="data row0 col4">0.67</td>
      <td id="T__row0_col6" class="data row0 col6">0.87</td>
      <td id="T__row0_col8" class="data row0 col8">0.84</td>
      <td id="T__row0_col9" class="data row0 col9">0.07</td>
      <td id="T__row0_col10" class="data row0 col10">0.11</td>
      <td id="T__row0_col11" class="data row0 col11">0.61</td>
      <td id="T__row0_col14" class="data row0 col14">UIP</td>
    </tr>
    <tr>
      <th id="T__level0_row1" class="row_heading level0 row1">M3GNet</th>
      <td id="T__row1_col0" class="data row1 col0">0.58</td>
      <td id="T__row1_col1" class="data row1 col1">2.66</td>
      <td id="T__row1_col2" class="data row1 col2">0.45</td>
      <td id="T__row1_col4" class="data row1 col4">0.79</td>
      <td id="T__row1_col6" class="data row1 col6">0.80</td>
      <td id="T__row1_col8" class="data row1 col8">0.80</td>
      <td id="T__row1_col9" class="data row1 col9">0.07</td>
      <td id="T__row1_col10" class="data row1 col10">0.12</td>
      <td id="T__row1_col11" class="data row1 col11">0.59</td>
      <td id="T__row1_col14" class="data row1 col14">UIP</td>
    </tr>
    <tr>
      <th id="T__level0_row2" class="row_heading level0 row2">MEGNet</th>
      <td id="T__row2_col0" class="data row2 col0">0.52</td>
      <td id="T__row2_col1" class="data row2 col1">2.70</td>
      <td id="T__row2_col2" class="data row2 col2">0.46</td>
      <td id="T__row2_col4" class="data row2 col4">0.59</td>
      <td id="T__row2_col6" class="data row2 col6">0.86</td>
      <td id="T__row2_col8" class="data row2 col8">0.81</td>
      <td id="T__row2_col9" class="data row2 col9">0.13</td>
      <td id="T__row2_col10" class="data row2 col10">0.20</td>
      <td id="T__row2_col11" class="data row2 col11">-0.27</td>
      <td id="T__row2_col14" class="data row2 col14">GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row3" class="row_heading level0 row3">CGCNN</th>
      <td id="T__row3_col0" class="data row3 col0">0.52</td>
      <td id="T__row3_col1" class="data row3 col1">2.62</td>
      <td id="T__row3_col2" class="data row3 col2">0.45</td>
      <td id="T__row3_col4" class="data row3 col4">0.60</td>
      <td id="T__row3_col6" class="data row3 col6">0.85</td>
      <td id="T__row3_col8" class="data row3 col8">0.81</td>
      <td id="T__row3_col9" class="data row3 col9">0.14</td>
      <td id="T__row3_col10" class="data row3 col10">0.23</td>
      <td id="T__row3_col11" class="data row3 col11">-0.61</td>
      <td id="T__row3_col14" class="data row3 col14">GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row4" class="row_heading level0 row4">CGCNN+P</th>
      <td id="T__row4_col0" class="data row4 col0">0.51</td>
      <td id="T__row4_col1" class="data row4 col1">2.38</td>
      <td id="T__row4_col2" class="data row4 col2">0.41</td>
      <td id="T__row4_col4" class="data row4 col4">0.69</td>
      <td id="T__row4_col6" class="data row4 col6">0.79</td>
      <td id="T__row4_col8" class="data row4 col8">0.78</td>
      <td id="T__row4_col9" class="data row4 col9">0.11</td>
      <td id="T__row4_col10" class="data row4 col10">0.18</td>
      <td id="T__row4_col11" class="data row4 col11">0.02</td>
      <td id="T__row4_col14" class="data row4 col14">GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row5" class="row_heading level0 row5">Wrenformer</th>
      <td id="T__row5_col0" class="data row5 col0">0.48</td>
      <td id="T__row5_col1" class="data row5 col1">2.13</td>
      <td id="T__row5_col2" class="data row5 col2">0.36</td>
      <td id="T__row5_col4" class="data row5 col4">0.71</td>
      <td id="T__row5_col6" class="data row5 col6">0.74</td>
      <td id="T__row5_col8" class="data row5 col8">0.74</td>
      <td id="T__row5_col9" class="data row5 col9">0.10</td>
      <td id="T__row5_col10" class="data row5 col10">0.18</td>
      <td id="T__row5_col11" class="data row5 col11">-0.04</td>
      <td id="T__row5_col14" class="data row5 col14">Transformer</td>
    </tr>
    <tr>
      <th id="T__level0_row6" class="row_heading level0 row6">BOWSR + MEGNet</th>
      <td id="T__row6_col0" class="data row6 col0">0.44</td>
      <td id="T__row6_col1" class="data row6 col1">1.90</td>
      <td id="T__row6_col2" class="data row6 col2">0.32</td>
      <td id="T__row6_col4" class="data row6 col4">0.74</td>
      <td id="T__row6_col6" class="data row6 col6">0.67</td>
      <td id="T__row6_col8" class="data row6 col8">0.68</td>
      <td id="T__row6_col9" class="data row6 col9">0.11</td>
      <td id="T__row6_col10" class="data row6 col10">0.16</td>
      <td id="T__row6_col11" class="data row6 col11">0.15</td>
      <td id="T__row6_col14" class="data row6 col14">BO+GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row7" class="row_heading level0 row7">Voronoi RF</th>
      <td id="T__row7_col0" class="data row7 col0">0.34</td>
      <td id="T__row7_col1" class="data row7 col1">1.51</td>
      <td id="T__row7_col2" class="data row7 col2">0.26</td>
      <td id="T__row7_col4" class="data row7 col4">0.52</td>
      <td id="T__row7_col6" class="data row7 col6">0.69</td>
      <td id="T__row7_col8" class="data row7 col8">0.66</td>
      <td id="T__row7_col9" class="data row7 col9">0.14</td>
      <td id="T__row7_col10" class="data row7 col10">0.21</td>
      <td id="T__row7_col11" class="data row7 col11">-0.32</td>
      <td id="T__row7_col14" class="data row7 col14">Fingerprint+RF</td>
    </tr>
    <tr>
      <th id="T__level0_row8" class="row_heading level0 row8">dummy</th>
      <td id="T__row8_col0" class="data row8 col0">0.19</td>
      <td id="T__row8_col1" class="data row8 col1">1.01</td>
      <td id="T__row8_col2" class="data row8 col2">0.17</td>
      <td id="T__row8_col4" class="data row8 col4">0.23</td>
      <td id="T__row8_col6" class="data row8 col6">0.77</td>
      <td id="T__row8_col8" class="data row8 col8">0.68</td>
      <td id="T__row8_col9" class="data row8 col9">0.12</td>
      <td id="T__row8_col10" class="data row8 col10">0.18</td>
      <td id="T__row8_col11" class="data row8 col11">0.00</td>
      <td id="T__row8_col14" class="data row8 col14">dummy clf</td>
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
  #T__row0_col6,
  #T__row0_col8,
  #T__row0_col9,
  #T__row0_col10,
  #T__row0_col11,
  #T__row1_col4 {
    background-color: #fde725;
    color: #000000;
  }
  #T__row0_col4 {
    background-color: #75d054;
    color: #000000;
  }
  #T__row1_col0 {
    background-color: #ece51b;
    color: #000000;
  }
  #T__row1_col1,
  #T__row4_col0 {
    background-color: #81d34d;
    color: #000000;
  }
  #T__row1_col2 {
    background-color: #7fd34e;
    color: #000000;
  }
  #T__row1_col6 {
    background-color: #32b67a;
    color: #f1f1f1;
  }
  #T__row1_col8,
  #T__row3_col1 {
    background-color: #73d056;
    color: #000000;
  }
  #T__row1_col9 {
    background-color: #f6e620;
    color: #000000;
  }
  #T__row1_col10 {
    background-color: #efe51c;
    color: #000000;
  }
  #T__row1_col11 {
    background-color: #f8e621;
    color: #000000;
  }
  #T__row2_col0 {
    background-color: #89d548;
    color: #000000;
  }
  #T__row2_col1 {
    background-color: #8bd646;
    color: #000000;
  }
  #T__row2_col2 {
    background-color: #90d743;
    color: #000000;
  }
  #T__row2_col4 {
    background-color: #2ab07f;
    color: #f1f1f1;
  }
  #T__row2_col6 {
    background-color: #d8e219;
    color: #000000;
  }
  #T__row2_col8,
  #T__row5_col4 {
    background-color: #a0da39;
    color: #000000;
  }
  #T__row2_col9 {
    background-color: #414287;
    color: #f1f1f1;
  }
  #T__row2_col10 {
    background-color: #3e4989;
    color: #f1f1f1;
  }
  #T__row2_col11 {
    background-color: #375b8d;
    color: #f1f1f1;
  }
  #T__row3_col0,
  #T__row4_col4 {
    background-color: #84d44b;
    color: #000000;
  }
  #T__row3_col2 {
    background-color: #7ad151;
    color: #000000;
  }
  #T__row3_col4 {
    background-color: #34b679;
    color: #f1f1f1;
  }
  #T__row3_col6 {
    background-color: #b0dd2f;
    color: #000000;
  }
  #T__row3_col8 {
    background-color: #86d549;
    color: #000000;
  }
  #T__row3_col9 {
    background-color: #482475;
    color: #f1f1f1;
  }
  #T__row3_col10,
  #T__row3_col11,
  #T__row6_col6,
  #T__row7_col8,
  #T__row7_col9,
  #T__row8_col0,
  #T__row8_col1,
  #T__row8_col2,
  #T__row8_col4 {
    background-color: #440154;
    color: #f1f1f1;
  }
  #T__row4_col1 {
    background-color: #38b977;
    color: #f1f1f1;
  }
  #T__row4_col2 {
    background-color: #3bbb75;
    color: #f1f1f1;
  }
  #T__row4_col6 {
    background-color: #25ac82;
    color: #f1f1f1;
  }
  #T__row4_col8 {
    background-color: #2fb47c;
    color: #f1f1f1;
  }
  #T__row4_col9 {
    background-color: #24878e;
    color: #f1f1f1;
  }
  #T__row4_col10 {
    background-color: #277f8e;
    color: #f1f1f1;
  }
  #T__row4_col11,
  #T__row5_col9 {
    background-color: #1f958b;
    color: #f1f1f1;
  }
  #T__row5_col0 {
    background-color: #54c568;
    color: #000000;
  }
  #T__row5_col1 {
    background-color: #1e9c89;
    color: #f1f1f1;
  }
  #T__row5_col2 {
    background-color: #1f9e89;
    color: #f1f1f1;
  }
  #T__row5_col6 {
    background-color: #2d718e;
    color: #f1f1f1;
  }
  #T__row5_col8 {
    background-color: #26818e;
    color: #f1f1f1;
  }
  #T__row5_col10 {
    background-color: #2b748e;
    color: #f1f1f1;
  }
  #T__row5_col11 {
    background-color: #23898e;
    color: #f1f1f1;
  }
  #T__row6_col0 {
    background-color: #29af7f;
    color: #f1f1f1;
  }
  #T__row6_col1 {
    background-color: #26828e;
    color: #f1f1f1;
  }
  #T__row6_col2,
  #T__row8_col10 {
    background-color: #287c8e;
    color: #f1f1f1;
  }
  #T__row6_col4 {
    background-color: #c8e020;
    color: #000000;
  }
  #T__row6_col8 {
    background-color: #482979;
    color: #f1f1f1;
  }
  #T__row6_col9 {
    background-color: #2a768e;
    color: #f1f1f1;
  }
  #T__row6_col10 {
    background-color: #1fa088;
    color: #f1f1f1;
  }
  #T__row6_col11 {
    background-color: #28ae80;
    color: #f1f1f1;
  }
  #T__row7_col0 {
    background-color: #2c738e;
    color: #f1f1f1;
  }
  #T__row7_col1,
  #T__row7_col11 {
    background-color: #3c508b;
    color: #f1f1f1;
  }
  #T__row7_col2 {
    background-color: #3a548c;
    color: #f1f1f1;
  }
  #T__row7_col4 {
    background-color: #20938c;
    color: #f1f1f1;
  }
  #T__row7_col6 {
    background-color: #482878;
    color: #f1f1f1;
  }
  #T__row7_col10 {
    background-color: #424086;
    color: #f1f1f1;
  }
  #T__row8_col6 {
    background-color: #21908d;
    color: #f1f1f1;
  }
  #T__row8_col8 {
    background-color: #482576;
    color: #f1f1f1;
  }
  #T__row8_col9 {
    background-color: #34618d;
    color: #f1f1f1;
  }
  #T__row8_col11 {
    background-color: #21918c;
    color: #f1f1f1;
  }
  #T_ :is(td, th):nth-last-child(3) {
    border-left: 1px dotted white;
  }
</style>
