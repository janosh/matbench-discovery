<script lang="ts">
  import { sortable } from 'svelte-zoo/actions'
</script>

<table use:sortable class="roomy" {...$$props} id="T_">
  <thead>
    <tr>
      <th class="index_name level0">Model</th>
      <th id="T__level0_col0" class="col_heading level0 col0">F1</th>
      <th id="T__level0_col1" class="col_heading level0 col1">DAF</th>
      <th id="T__level0_col2" class="col_heading level0 col2">Precision</th>
      <th id="T__level0_col3" class="col_heading level0 col3">Accuracy</th>
      <th id="T__level0_col4" class="col_heading level0 col4">MAE</th>
      <th id="T__level0_col5" class="col_heading level0 col5">RMSE</th>
      <th id="T__level0_col6" class="col_heading level0 col6">R<sup>2</sup></th>
      <th id="T__level0_col7" class="col_heading level0 col7">Model Class</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T__level0_row0" class="row_heading level0 row0">CHGNet</th>
      <td id="T__row0_col0" class="data row0 col0">0.91</td>
      <td id="T__row0_col1" class="data row0 col1">4.97</td>
      <td id="T__row0_col2" class="data row0 col2">0.83</td>
      <td id="T__row0_col3" class="data row0 col3">0.83</td>
      <td id="T__row0_col4" class="data row0 col4">0.07</td>
      <td id="T__row0_col5" class="data row0 col5">0.10</td>
      <td id="T__row0_col6" class="data row0 col6">0.82</td>
      <td id="T__row0_col7" class="data row0 col7">UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row1" class="row_heading level0 row1">M3GNet</th>
      <td id="T__row1_col0" class="data row1 col0">0.89</td>
      <td id="T__row1_col1" class="data row1 col1">4.77</td>
      <td id="T__row1_col2" class="data row1 col2">0.80</td>
      <td id="T__row1_col3" class="data row1 col3">0.80</td>
      <td id="T__row1_col4" class="data row1 col4">0.10</td>
      <td id="T__row1_col5" class="data row1 col5">0.16</td>
      <td id="T__row1_col6" class="data row1 col6">0.61</td>
      <td id="T__row1_col7" class="data row1 col7">UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row2" class="row_heading level0 row2">CGCNN</th>
      <td id="T__row2_col0" class="data row2 col0">0.76</td>
      <td id="T__row2_col1" class="data row2 col1">3.66</td>
      <td id="T__row2_col2" class="data row2 col2">0.61</td>
      <td id="T__row2_col3" class="data row2 col3">0.61</td>
      <td id="T__row2_col4" class="data row2 col4">0.17</td>
      <td id="T__row2_col5" class="data row2 col5">0.24</td>
      <td id="T__row2_col6" class="data row2 col6">0.26</td>
      <td id="T__row2_col7" class="data row2 col7">GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row3" class="row_heading level0 row3">CGCNN+P</th>
      <td id="T__row3_col0" class="data row3 col0">0.74</td>
      <td id="T__row3_col1" class="data row3 col1">3.51</td>
      <td id="T__row3_col2" class="data row3 col2">0.58</td>
      <td id="T__row3_col3" class="data row3 col3">0.58</td>
      <td id="T__row3_col4" class="data row3 col4">0.21</td>
      <td id="T__row3_col5" class="data row3 col5">0.29</td>
      <td id="T__row3_col6" class="data row3 col6">0.02</td>
      <td id="T__row3_col7" class="data row3 col7">GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row4" class="row_heading level0 row4">Wrenformer</th>
      <td id="T__row4_col0" class="data row4 col0">0.74</td>
      <td id="T__row4_col1" class="data row4 col1">3.51</td>
      <td id="T__row4_col2" class="data row4 col2">0.58</td>
      <td id="T__row4_col3" class="data row4 col3">0.58</td>
      <td id="T__row4_col4" class="data row4 col4">0.19</td>
      <td id="T__row4_col5" class="data row4 col5">0.25</td>
      <td id="T__row4_col6" class="data row4 col6">0.21</td>
      <td id="T__row4_col7" class="data row4 col7">Transformer</td>
    </tr>
    <tr>
      <th id="T__level0_row5" class="row_heading level0 row5">BOWSR + MEGNet</th>
      <td id="T__row5_col0" class="data row5 col0">0.69</td>
      <td id="T__row5_col1" class="data row5 col1">3.14</td>
      <td id="T__row5_col2" class="data row5 col2">0.52</td>
      <td id="T__row5_col3" class="data row5 col3">0.52</td>
      <td id="T__row5_col4" class="data row5 col4">0.27</td>
      <td id="T__row5_col5" class="data row5 col5">0.33</td>
      <td id="T__row5_col6" class="data row5 col6">-1.10</td>
      <td id="T__row5_col7" class="data row5 col7">BO+GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row6" class="row_heading level0 row6">MEGNet</th>
      <td id="T__row6_col0" class="data row6 col0">0.65</td>
      <td id="T__row6_col1" class="data row6 col1">2.86</td>
      <td id="T__row6_col2" class="data row6 col2">0.48</td>
      <td id="T__row6_col3" class="data row6 col3">0.48</td>
      <td id="T__row6_col4" class="data row6 col4">0.33</td>
      <td id="T__row6_col5" class="data row6 col5">0.36</td>
      <td id="T__row6_col6" class="data row6 col6">-0.82</td>
      <td id="T__row6_col7" class="data row6 col7">GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row7" class="row_heading level0 row7">Voronoi RF</th>
      <td id="T__row7_col0" class="data row7 col0">0.58</td>
      <td id="T__row7_col1" class="data row7 col1">2.47</td>
      <td id="T__row7_col2" class="data row7 col2">0.41</td>
      <td id="T__row7_col3" class="data row7 col3">0.41</td>
      <td id="T__row7_col4" class="data row7 col4">0.36</td>
      <td id="T__row7_col5" class="data row7 col5">0.43</td>
      <td id="T__row7_col6" class="data row7 col6">-0.83</td>
      <td id="T__row7_col7" class="data row7 col7">Fingerprint+RF</td>
    </tr>
    <tr>
      <th id="T__level0_row8" class="row_heading level0 row8">Dummy</th>
      <td id="T__row8_col0" class="data row8 col0">0.19</td>
      <td id="T__row8_col1" class="data row8 col1">1.00</td>
      <td id="T__row8_col2" class="data row8 col2">0.17</td>
      <td id="T__row8_col3" class="data row8 col3">0.68</td>
      <td id="T__row8_col4" class="data row8 col4">0.12</td>
      <td id="T__row8_col5" class="data row8 col5">0.18</td>
      <td id="T__row8_col6" class="data row8 col6">0.00</td>
      <td id="T__row8_col7" class="data row8 col7">scikit-learn</td>
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
  #T__row0_col4,
  #T__row0_col5,
  #T__row0_col6 {
    background-color: #fde725;
    color: #000000;
  }
  #T__row1_col0 {
    background-color: #ece51b;
    color: #000000;
  }
  #T__row1_col1,
  #T__row1_col2 {
    background-color: #dfe318;
    color: #000000;
  }
  #T__row1_col3 {
    background-color: #cae11f;
    color: #000000;
  }
  #T__row1_col4 {
    background-color: #bddf26;
    color: #000000;
  }
  #T__row1_col5 {
    background-color: #89d548;
    color: #000000;
  }
  #T__row1_col6 {
    background-color: #b8de29;
    color: #000000;
  }
  #T__row2_col0 {
    background-color: #75d054;
    color: #000000;
  }
  #T__row2_col1,
  #T__row2_col2 {
    background-color: #37b878;
    color: #f1f1f1;
  }
  #T__row2_col3 {
    background-color: #238a8d;
    color: #f1f1f1;
  }
  #T__row2_col4 {
    background-color: #2fb47c;
    color: #f1f1f1;
  }
  #T__row2_col5 {
    background-color: #21a585;
    color: #f1f1f1;
  }
  #T__row2_col6 {
    background-color: #4ac16d;
    color: #000000;
  }
  #T__row3_col0,
  #T__row4_col0 {
    background-color: #65cb5e;
    color: #000000;
  }
  #T__row3_col1,
  #T__row4_col1,
  #T__row6_col0 {
    background-color: #2ab07f;
    color: #f1f1f1;
  }
  #T__row3_col2,
  #T__row4_col2 {
    background-color: #29af7f;
    color: #f1f1f1;
  }
  #T__row3_col3,
  #T__row4_col3 {
    background-color: #287c8e;
    color: #f1f1f1;
  }
  #T__row3_col4 {
    background-color: #20928c;
    color: #f1f1f1;
  }
  #T__row3_col5 {
    background-color: #277f8e;
    color: #f1f1f1;
  }
  #T__row3_col6 {
    background-color: #20a486;
    color: #f1f1f1;
  }
  #T__row4_col4 {
    background-color: #21a685;
    color: #f1f1f1;
  }
  #T__row4_col5,
  #T__row5_col1 {
    background-color: #1f9a8a;
    color: #f1f1f1;
  }
  #T__row4_col6 {
    background-color: #3bbb75;
    color: #f1f1f1;
  }
  #T__row5_col0 {
    background-color: #40bd72;
    color: #f1f1f1;
  }
  #T__row5_col2 {
    background-color: #1f998a;
    color: #f1f1f1;
  }
  #T__row5_col3 {
    background-color: #39568c;
    color: #f1f1f1;
  }
  #T__row5_col4 {
    background-color: #31678e;
    color: #f1f1f1;
  }
  #T__row5_col5 {
    background-color: #33638d;
    color: #f1f1f1;
  }
  #T__row5_col6,
  #T__row7_col3,
  #T__row7_col4,
  #T__row7_col5,
  #T__row8_col0,
  #T__row8_col1,
  #T__row8_col2 {
    background-color: #440154;
    color: #f1f1f1;
  }
  #T__row6_col1,
  #T__row6_col2 {
    background-color: #23888e;
    color: #f1f1f1;
  }
  #T__row6_col3 {
    background-color: #453581;
    color: #f1f1f1;
  }
  #T__row6_col4 {
    background-color: #482576;
    color: #f1f1f1;
  }
  #T__row6_col5 {
    background-color: #3e4a89;
    color: #f1f1f1;
  }
  #T__row6_col6 {
    background-color: #46337f;
    color: #f1f1f1;
  }
  #T__row7_col0 {
    background-color: #1e9b8a;
    color: #f1f1f1;
  }
  #T__row7_col1 {
    background-color: #2c718e;
    color: #f1f1f1;
  }
  #T__row7_col2 {
    background-color: #2d718e;
    color: #f1f1f1;
  }
  #T__row7_col6 {
    background-color: #46327e;
    color: #f1f1f1;
  }
  #T__row8_col3 {
    background-color: #2eb37c;
    color: #f1f1f1;
  }
  #T__row8_col4 {
    background-color: #90d743;
    color: #000000;
  }
  #T__row8_col5 {
    background-color: #67cc5c;
    color: #000000;
  }
  #T__row8_col6 {
    background-color: #1fa287;
    color: #f1f1f1;
  }
  #T_ :is(td, th):nth-last-child(4) {
    border-left: 1px dotted white;
  }
</style>
