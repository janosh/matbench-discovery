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
#T_ th {
  border: 1px solid;
  border-width: 1px 0;
  text-align: left;
}
#T__row0_col0, #T__row0_col1, #T__row0_col2, #T__row0_col3, #T__row0_col4, #T__row0_col5, #T__row0_col6 {
  background-color: #fde725;
  color: #000000;
}
#T__row1_col0 {
  background-color: #d2e21b;
  color: #000000;
}
#T__row1_col1, #T__row1_col2, #T__row1_col3 {
  background-color: #cae11f;
  color: #000000;
}
#T__row1_col4 {
  background-color: #c0df25;
  color: #000000;
}
#T__row1_col5 {
  background-color: #7ad151;
  color: #000000;
}
#T__row1_col6 {
  background-color: #bade28;
  color: #000000;
}
#T__row2_col0, #T__row2_col4 {
  background-color: #470d60;
  color: #f1f1f1;
}
#T__row2_col1, #T__row2_col6 {
  background-color: #460a5d;
  color: #f1f1f1;
}
#T__row2_col2, #T__row2_col3 {
  background-color: #460b5e;
  color: #f1f1f1;
}
#T__row2_col5 {
  background-color: #471164;
  color: #f1f1f1;
}
#T__row3_col0, #T__row3_col1, #T__row3_col2, #T__row3_col3, #T__row4_col2, #T__row4_col3, #T__row5_col6 {
  background-color: #450457;
  color: #f1f1f1;
}
#T__row3_col4 {
  background-color: #482173;
  color: #f1f1f1;
}
#T__row3_col5 {
  background-color: #482071;
  color: #f1f1f1;
}
#T__row3_col6 {
  background-color: #46337f;
  color: #f1f1f1;
}
#T__row4_col0, #T__row4_col1 {
  background-color: #440256;
  color: #f1f1f1;
}
#T__row4_col4, #T__row4_col5, #T__row4_col6, #T__row5_col0, #T__row5_col1, #T__row5_col2, #T__row5_col3, #T__row5_col4, #T__row5_col5 {
  background-color: #440154;
  color: #f1f1f1;
}
#T_ :is(td, th):nth-child(6) { border-left: 1px dotted white; }
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

    <table use:sortable class='roomy' {...$$props}
     id="T_">
  <thead>
    <tr>
      <th class="index_name level0" >Model</th>
      <th id="T__level0_col0" class="col_heading level0 col0" >F1</th>
      <th id="T__level0_col1" class="col_heading level0 col1" >DAF</th>
      <th id="T__level0_col2" class="col_heading level0 col2" >Precision</th>
      <th id="T__level0_col3" class="col_heading level0 col3" >Accuracy</th>
      <th id="T__level0_col4" class="col_heading level0 col4" >MAE</th>
      <th id="T__level0_col5" class="col_heading level0 col5" >RMSE</th>
      <th id="T__level0_col6" class="col_heading level0 col6" >R<sup>2</sup></th>
      <th id="T__level0_col7" class="col_heading level0 col7" >training size</th>
      <th id="T__level0_col8" class="col_heading level0 col8" >Model Class</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T__level0_row0" class="row_heading level0 row0" >CHGNet</th>
      <td id="T__row0_col0" class="data row0 col0" >0.91</td>
      <td id="T__row0_col1" class="data row0 col1" >4.96</td>
      <td id="T__row0_col2" class="data row0 col2" >0.83</td>
      <td id="T__row0_col3" class="data row0 col3" >0.83</td>
      <td id="T__row0_col4" class="data row0 col4" >0.07</td>
      <td id="T__row0_col5" class="data row0 col5" >0.10</td>
      <td id="T__row0_col6" class="data row0 col6" >0.82</td>
      <td id="T__row0_col7" class="data row0 col7" >1,580,395</td>
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
      <td id="T__row1_col7" class="data row1 col7" >188,349</td>
      <td id="T__row1_col8" class="data row1 col8" >UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row2" class="row_heading level0 row2" >MEGNet RS2RE</th>
      <td id="T__row2_col0" class="data row2 col0" >0.65</td>
      <td id="T__row2_col1" class="data row2 col1" >2.88</td>
      <td id="T__row2_col2" class="data row2 col2" >0.48</td>
      <td id="T__row2_col3" class="data row2 col3" >0.48</td>
      <td id="T__row2_col4" class="data row2 col4" >0.35</td>
      <td id="T__row2_col5" class="data row2 col5" >0.37</td>
      <td id="T__row2_col6" class="data row2 col6" >-1.06</td>
      <td id="T__row2_col7" class="data row2 col7" >69,239</td>
      <td id="T__row2_col8" class="data row2 col8" >GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row3" class="row_heading level0 row3" >MEGNet</th>
      <td id="T__row3_col0" class="data row3 col0" >0.64</td>
      <td id="T__row3_col1" class="data row3 col1" >2.84</td>
      <td id="T__row3_col2" class="data row3 col2" >0.47</td>
      <td id="T__row3_col3" class="data row3 col3" >0.47</td>
      <td id="T__row3_col4" class="data row3 col4" >0.33</td>
      <td id="T__row3_col5" class="data row3 col5" >0.36</td>
      <td id="T__row3_col6" class="data row3 col6" >-0.82</td>
      <td id="T__row3_col7" class="data row3 col7" >69,239</td>
      <td id="T__row3_col8" class="data row3 col8" >GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row4" class="row_heading level0 row4" >CHGNet→MEGNet</th>
      <td id="T__row4_col0" class="data row4 col0" >0.64</td>
      <td id="T__row4_col1" class="data row4 col1" >2.83</td>
      <td id="T__row4_col2" class="data row4 col2" >0.47</td>
      <td id="T__row4_col3" class="data row4 col3" >0.47</td>
      <td id="T__row4_col4" class="data row4 col4" >0.36</td>
      <td id="T__row4_col5" class="data row4 col5" >0.38</td>
      <td id="T__row4_col6" class="data row4 col6" >-1.11</td>
      <td id="T__row4_col7" class="data row4 col7" >1,580,395</td>
      <td id="T__row4_col8" class="data row4 col8" >UIP-GNN</td>
    </tr>
    <tr>
      <th id="T__level0_row5" class="row_heading level0 row5" >M3GNet→MEGNet</th>
      <td id="T__row5_col0" class="data row5 col0" >0.64</td>
      <td id="T__row5_col1" class="data row5 col1" >2.82</td>
      <td id="T__row5_col2" class="data row5 col2" >0.47</td>
      <td id="T__row5_col3" class="data row5 col3" >0.47</td>
      <td id="T__row5_col4" class="data row5 col4" >0.36</td>
      <td id="T__row5_col5" class="data row5 col5" >0.38</td>
      <td id="T__row5_col6" class="data row5 col6" >-1.09</td>
      <td id="T__row5_col7" class="data row5 col7" >188,349</td>
      <td id="T__row5_col8" class="data row5 col8" >UIP-GNN</td>
    </tr>
  </tbody>
</table>
