<style type="text/css">
#T_  {
  font-family: sans-serif;
  border-collapse: collapse;
}
#T_ td {
  border: none;
  padding: 4px 6px;
  white-space: nowrap;
}
#T_  th {
  border: none;
  padding: 4px 6px;
  white-space: nowrap;
}
#T_ th.col_heading {
  border: 1px solid;
  border-width: 1px 0;
  text-align: left;
}
#T_ th.row_heading {
  font-weight: normal;
  padding: 3pt;
}
#T__row0_col0, #T__row0_col1, #T__row0_col2, #T__row0_col3, #T__row0_col6, #T__row0_col7, #T__row0_col8, #T__row1_col4, #T__row4_col5 {
  background-color: #fde725;
  color: #000000;
}
#T__row0_col4 {
  background-color: #77d153;
  color: #000000;
}
#T__row0_col5 {
  background-color: #eae51a;
  color: #000000;
}
#T__row1_col0 {
  background-color: #1fa187;
  color: #f1f1f1;
}
#T__row1_col1 {
  background-color: #433d84;
  color: #f1f1f1;
}
#T__row1_col2 {
  background-color: #404588;
  color: #f1f1f1;
}
#T__row1_col3 {
  background-color: #3f4889;
  color: #f1f1f1;
}
#T__row1_col5 {
  background-color: #481d6f;
  color: #f1f1f1;
}
#T__row1_col6, #T__row2_col4 {
  background-color: #89d548;
  color: #000000;
}
#T__row1_col7 {
  background-color: #9dd93b;
  color: #000000;
}
#T__row1_col8 {
  background-color: #b5de2b;
  color: #000000;
}
#T__row2_col0 {
  background-color: #31678e;
  color: #f1f1f1;
}
#T__row2_col1, #T__row2_col2, #T__row2_col5, #T__row3_col2, #T__row3_col3, #T__row4_col0, #T__row4_col4, #T__row4_col6, #T__row4_col7, #T__row4_col8 {
  background-color: #440154;
  color: #f1f1f1;
}
#T__row2_col3 {
  background-color: #46085c;
  color: #f1f1f1;
}
#T__row2_col6, #T__row3_col7 {
  background-color: #20a486;
  color: #f1f1f1;
}
#T__row2_col7 {
  background-color: #21a685;
  color: #f1f1f1;
}
#T__row2_col8 {
  background-color: #38b977;
  color: #f1f1f1;
}
#T__row3_col0 {
  background-color: #34608d;
  color: #f1f1f1;
}
#T__row3_col1 {
  background-color: #460b5e;
  color: #f1f1f1;
}
#T__row3_col4 {
  background-color: #63cb5f;
  color: #000000;
}
#T__row3_col5 {
  background-color: #470d60;
  color: #f1f1f1;
}
#T__row3_col6 {
  background-color: #1e9d89;
  color: #f1f1f1;
}
#T__row3_col8 {
  background-color: #31b57b;
  color: #f1f1f1;
}
#T__row4_col1 {
  background-color: #33628d;
  color: #f1f1f1;
}
#T__row4_col2 {
  background-color: #31668e;
  color: #f1f1f1;
}
#T__row4_col3 {
  background-color: #20928c;
  color: #f1f1f1;
}
#T_ :is(td, th):nth-child(8) { border-left: 1px solid white; }
    table {
        scrollbar-width: none;  /* Firefox */
    }
    table::-webkit-scrollbar {
        display: none;  /* Safari and Chrome */
    }
</style>
<script lang="ts">
      import { sortable } from 'svelte-zoo/actions'
    </script>

    <table class='roomy' use:sortable {...$$props}
     id="T_">
  <thead>
    <tr>
      <th class="index_name level0" >Model</th>
      <th id="T__level0_col0" class="col_heading level0 col0" >F1 ↑</th>
      <th id="T__level0_col1" class="col_heading level0 col1" >DAF ↑</th>
      <th id="T__level0_col2" class="col_heading level0 col2" >Prec ↑</th>
      <th id="T__level0_col3" class="col_heading level0 col3" >Acc ↑</th>
      <th id="T__level0_col4" class="col_heading level0 col4" >TPR ↑</th>
      <th id="T__level0_col5" class="col_heading level0 col5" >TNR ↑</th>
      <th id="T__level0_col6" class="col_heading level0 col6" >MAE ↓</th>
      <th id="T__level0_col7" class="col_heading level0 col7" >RMSE ↓</th>
      <th id="T__level0_col8" class="col_heading level0 col8" >R<sup>2</sup> ↑</th>
      <th id="T__level0_col9" class="col_heading level0 col9" >Training Size</th>
      <th id="T__level0_col10" class="col_heading level0 col10" >Model Type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T__level0_row0" class="row_heading level0 row0" >CHGNet</th>
      <td id="T__row0_col0" class="data row0 col0" >0.61</td>
      <td id="T__row0_col1" class="data row0 col1" >3.32</td>
      <td id="T__row0_col2" class="data row0 col2" >0.51</td>
      <td id="T__row0_col3" class="data row0 col3" >0.85</td>
      <td id="T__row0_col4" class="data row0 col4" >0.76</td>
      <td id="T__row0_col5" class="data row0 col5" >0.87</td>
      <td id="T__row0_col6" class="data row0 col6" >0.06</td>
      <td id="T__row0_col7" class="data row0 col7" >0.10</td>
      <td id="T__row0_col8" class="data row0 col8" >0.69</td>
      <td id="T__row0_col9" class="data row0 col9" >1.6M <small>(145.9K)</small></td>
      <td id="T__row0_col10" class="data row0 col10" >UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row1" class="row_heading level0 row1" >M3GNet</th>
      <td id="T__row1_col0" class="data row1 col0" >0.57</td>
      <td id="T__row1_col1" class="data row1 col1" >2.84</td>
      <td id="T__row1_col2" class="data row1 col2" >0.44</td>
      <td id="T__row1_col3" class="data row1 col3" >0.81</td>
      <td id="T__row1_col4" class="data row1 col4" >0.80</td>
      <td id="T__row1_col5" class="data row1 col5" >0.81</td>
      <td id="T__row1_col6" class="data row1 col6" >0.07</td>
      <td id="T__row1_col7" class="data row1 col7" >0.12</td>
      <td id="T__row1_col8" class="data row1 col8" >0.58</td>
      <td id="T__row1_col9" class="data row1 col9" >188.3K <small>(62.8K)</small></td>
      <td id="T__row1_col10" class="data row1 col10" >UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row2" class="row_heading level0 row2" >M3GNet→MEGNet</th>
      <td id="T__row2_col0" class="data row2 col0" >0.54</td>
      <td id="T__row2_col1" class="data row2 col1" >2.74</td>
      <td id="T__row2_col2" class="data row2 col2" >0.42</td>
      <td id="T__row2_col3" class="data row2 col3" >0.80</td>
      <td id="T__row2_col4" class="data row2 col4" >0.76</td>
      <td id="T__row2_col5" class="data row2 col5" >0.81</td>
      <td id="T__row2_col6" class="data row2 col6" >0.09</td>
      <td id="T__row2_col7" class="data row2 col7" >0.14</td>
      <td id="T__row2_col8" class="data row2 col8" >0.38</td>
      <td id="T__row2_col9" class="data row2 col9" >188.3K <small>(62.8K)</small></td>
      <td id="T__row2_col10" class="data row2 col10" >UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row3" class="row_heading level0 row3" >CHGNet→MEGNet</th>
      <td id="T__row3_col0" class="data row3 col0" >0.54</td>
      <td id="T__row3_col1" class="data row3 col1" >2.76</td>
      <td id="T__row3_col2" class="data row3 col2" >0.42</td>
      <td id="T__row3_col3" class="data row3 col3" >0.80</td>
      <td id="T__row3_col4" class="data row3 col4" >0.75</td>
      <td id="T__row3_col5" class="data row3 col5" >0.81</td>
      <td id="T__row3_col6" class="data row3 col6" >0.09</td>
      <td id="T__row3_col7" class="data row3 col7" >0.15</td>
      <td id="T__row3_col8" class="data row3 col8" >0.36</td>
      <td id="T__row3_col9" class="data row3 col9" >1.6M <small>(145.9K)</small></td>
      <td id="T__row3_col10" class="data row3 col10" >UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row4" class="row_heading level0 row4" >MEGNet</th>
      <td id="T__row4_col0" class="data row4 col0" >0.51</td>
      <td id="T__row4_col1" class="data row4 col1" >2.92</td>
      <td id="T__row4_col2" class="data row4 col2" >0.45</td>
      <td id="T__row4_col3" class="data row4 col3" >0.83</td>
      <td id="T__row4_col4" class="data row4 col4" >0.58</td>
      <td id="T__row4_col5" class="data row4 col5" >0.87</td>
      <td id="T__row4_col6" class="data row4 col6" >0.13</td>
      <td id="T__row4_col7" class="data row4 col7" >0.21</td>
      <td id="T__row4_col8" class="data row4 col8" >-0.25</td>
      <td id="T__row4_col9" class="data row4 col9" >133.4K</td>
      <td id="T__row4_col10" class="data row4 col10" >GNN</td>
    </tr>
  </tbody>
</table>
