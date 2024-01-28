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
#T__row0_col0, #T__row0_col1, #T__row0_col2, #T__row0_col3, #T__row0_col4, #T__row0_col5, #T__row0_col6 {
  background-color: #fde725;
  color: #000000;
}
#T__row1_col0 {
  background-color: #90d743;
  color: #000000;
}
#T__row1_col1 {
  background-color: #7ad151;
  color: #000000;
}
#T__row1_col2, #T__row1_col3 {
  background-color: #77d153;
  color: #000000;
}
#T__row1_col4 {
  background-color: #a8db34;
  color: #000000;
}
#T__row1_col5 {
  background-color: #81d34d;
  color: #000000;
}
#T__row1_col6 {
  background-color: #bddf26;
  color: #000000;
}
#T__row2_col0, #T__row2_col2, #T__row2_col3, #T__row4_col6 {
  background-color: #450457;
  color: #f1f1f1;
}
#T__row2_col1, #T__row3_col0, #T__row3_col1, #T__row3_col2, #T__row3_col3 {
  background-color: #440256;
  color: #f1f1f1;
}
#T__row2_col4, #T__row2_col5 {
  background-color: #482071;
  color: #f1f1f1;
}
#T__row2_col6 {
  background-color: #46337f;
  color: #f1f1f1;
}
#T__row3_col4, #T__row3_col5, #T__row3_col6, #T__row4_col0, #T__row4_col1, #T__row4_col2, #T__row4_col3, #T__row4_col4, #T__row4_col5 {
  background-color: #440154;
  color: #f1f1f1;
}
#T_ :is(td, th):nth-child(6) { border-left: 1px solid white; }
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
      <th id="T__level0_col4" class="col_heading level0 col4" >MAE ↓</th>
      <th id="T__level0_col5" class="col_heading level0 col5" >RMSE ↓</th>
      <th id="T__level0_col6" class="col_heading level0 col6" >R<sup>2</sup> ↑</th>
      <th id="T__level0_col7" class="col_heading level0 col7" >Training Size</th>
      <th id="T__level0_col8" class="col_heading level0 col8" >Model Type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T__level0_row0" class="row_heading level0 row0" >CHGNet</th>
      <td id="T__row0_col0" class="data row0 col0" >0.94</td>
      <td id="T__row0_col1" class="data row0 col1" >5.29</td>
      <td id="T__row0_col2" class="data row0 col2" >0.88</td>
      <td id="T__row0_col3" class="data row0 col3" >0.88</td>
      <td id="T__row0_col4" class="data row0 col4" >0.06</td>
      <td id="T__row0_col5" class="data row0 col5" >0.11</td>
      <td id="T__row0_col6" class="data row0 col6" >0.81</td>
      <td id="T__row0_col7" class="data row0 col7" >1.6M <small>(145.9K)</small></td>
      <td id="T__row0_col8" class="data row0 col8" >UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row1" class="row_heading level0 row1" >M3GNet</th>
      <td id="T__row1_col0" class="data row1 col0" >0.89</td>
      <td id="T__row1_col1" class="data row1 col1" >4.79</td>
      <td id="T__row1_col2" class="data row1 col2" >0.80</td>
      <td id="T__row1_col3" class="data row1 col3" >0.80</td>
      <td id="T__row1_col4" class="data row1 col4" >0.10</td>
      <td id="T__row1_col5" class="data row1 col5" >0.16</td>
      <td id="T__row1_col6" class="data row1 col6" >0.62</td>
      <td id="T__row1_col7" class="data row1 col7" >188.3K <small>(62.8K)</small></td>
      <td id="T__row1_col8" class="data row1 col8" >UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row2" class="row_heading level0 row2" >MEGNet</th>
      <td id="T__row2_col0" class="data row2 col0" >0.64</td>
      <td id="T__row2_col1" class="data row2 col1" >2.84</td>
      <td id="T__row2_col2" class="data row2 col2" >0.47</td>
      <td id="T__row2_col3" class="data row2 col3" >0.47</td>
      <td id="T__row2_col4" class="data row2 col4" >0.33</td>
      <td id="T__row2_col5" class="data row2 col5" >0.36</td>
      <td id="T__row2_col6" class="data row2 col6" >-0.82</td>
      <td id="T__row2_col7" class="data row2 col7" >133.4K</td>
      <td id="T__row2_col8" class="data row2 col8" >GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row3" class="row_heading level0 row3" >CHGNet→MEGNet</th>
      <td id="T__row3_col0" class="data row3 col0" >0.64</td>
      <td id="T__row3_col1" class="data row3 col1" >2.83</td>
      <td id="T__row3_col2" class="data row3 col2" >0.47</td>
      <td id="T__row3_col3" class="data row3 col3" >0.47</td>
      <td id="T__row3_col4" class="data row3 col4" >0.36</td>
      <td id="T__row3_col5" class="data row3 col5" >0.38</td>
      <td id="T__row3_col6" class="data row3 col6" >-1.11</td>
      <td id="T__row3_col7" class="data row3 col7" >1.6M <small>(145.9K)</small></td>
      <td id="T__row3_col8" class="data row3 col8" >UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row4" class="row_heading level0 row4" >M3GNet→MEGNet</th>
      <td id="T__row4_col0" class="data row4 col0" >0.64</td>
      <td id="T__row4_col1" class="data row4 col1" >2.82</td>
      <td id="T__row4_col2" class="data row4 col2" >0.47</td>
      <td id="T__row4_col3" class="data row4 col3" >0.47</td>
      <td id="T__row4_col4" class="data row4 col4" >0.36</td>
      <td id="T__row4_col5" class="data row4 col5" >0.38</td>
      <td id="T__row4_col6" class="data row4 col6" >-1.09</td>
      <td id="T__row4_col7" class="data row4 col7" >188.3K <small>(62.8K)</small></td>
      <td id="T__row4_col8" class="data row4 col8" >UIP-GNN</td>
    </tr>
  </tbody>
</table>
