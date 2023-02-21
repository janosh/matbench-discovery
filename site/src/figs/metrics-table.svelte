<script>
  import { sortable } from 'svelte-zoo/actions'
</script>

<table use:sortable {...$$props} id="T_">
  <thead>
    <tr>
      <th class="index_name level0">model</th>
      <th id="T__level0_col0" class="col_heading level0 col0">F1</th>
      <th id="T__level0_col1" class="col_heading level0 col1">R²</th>
      <th id="T__level0_col2" class="col_heading level0 col2">DAF</th>
      <th id="T__level0_col3" class="col_heading level0 col3">Precision</th>
      <th id="T__level0_col4" class="col_heading level0 col4">Recall</th>
      <th id="T__level0_col5" class="col_heading level0 col5">Accuracy</th>
      <th id="T__level0_col6" class="col_heading level0 col6">TPR</th>
      <th id="T__level0_col7" class="col_heading level0 col7">FPR</th>
      <th id="T__level0_col8" class="col_heading level0 col8">TNR</th>
      <th id="T__level0_col9" class="col_heading level0 col9">FNR</th>
      <th id="T__level0_col10" class="col_heading level0 col10">MAE</th>
      <th id="T__level0_col11" class="col_heading level0 col11">RMSE</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T__level0_row0" class="row_heading level0 row0">Voronoi Random Forest</th>
      <td id="T__row0_col0" class="data row0 col0">0.34</td>
      <td id="T__row0_col1" class="data row0 col1">-0.32</td>
      <td id="T__row0_col2" class="data row0 col2">1.51</td>
      <td id="T__row0_col3" class="data row0 col3">0.26</td>
      <td id="T__row0_col4" class="data row0 col4">0.52</td>
      <td id="T__row0_col5" class="data row0 col5">0.66</td>
      <td id="T__row0_col6" class="data row0 col6">0.52</td>
      <td id="T__row0_col7" class="data row0 col7">0.31</td>
      <td id="T__row0_col8" class="data row0 col8">0.69</td>
      <td id="T__row0_col9" class="data row0 col9">0.48</td>
      <td id="T__row0_col10" class="data row0 col10">0.14</td>
      <td id="T__row0_col11" class="data row0 col11">0.21</td>
    </tr>
    <tr>
      <th id="T__level0_row1" class="row_heading level0 row1">BOWSR + MEGNet</th>
      <td id="T__row1_col0" class="data row1 col0">0.44</td>
      <td id="T__row1_col1" class="data row1 col1">0.15</td>
      <td id="T__row1_col2" class="data row1 col2">1.90</td>
      <td id="T__row1_col3" class="data row1 col3">0.32</td>
      <td id="T__row1_col4" class="data row1 col4">0.74</td>
      <td id="T__row1_col5" class="data row1 col5">0.68</td>
      <td id="T__row1_col6" class="data row1 col6">0.74</td>
      <td id="T__row1_col7" class="data row1 col7">0.33</td>
      <td id="T__row1_col8" class="data row1 col8">0.67</td>
      <td id="T__row1_col9" class="data row1 col9">0.26</td>
      <td id="T__row1_col10" class="data row1 col10">0.11</td>
      <td id="T__row1_col11" class="data row1 col11">0.16</td>
    </tr>
    <tr>
      <th id="T__level0_row2" class="row_heading level0 row2">Wrenformer</th>
      <td id="T__row2_col0" class="data row2 col0">0.48</td>
      <td id="T__row2_col1" class="data row2 col1">-0.04</td>
      <td id="T__row2_col2" class="data row2 col2">2.13</td>
      <td id="T__row2_col3" class="data row2 col3">0.36</td>
      <td id="T__row2_col4" class="data row2 col4">0.71</td>
      <td id="T__row2_col5" class="data row2 col5">0.74</td>
      <td id="T__row2_col6" class="data row2 col6">0.71</td>
      <td id="T__row2_col7" class="data row2 col7">0.26</td>
      <td id="T__row2_col8" class="data row2 col8">0.74</td>
      <td id="T__row2_col9" class="data row2 col9">0.29</td>
      <td id="T__row2_col10" class="data row2 col10">0.10</td>
      <td id="T__row2_col11" class="data row2 col11">0.18</td>
    </tr>
    <tr>
      <th id="T__level0_row3" class="row_heading level0 row3">MEGNet</th>
      <td id="T__row3_col0" class="data row3 col0">0.49</td>
      <td id="T__row3_col1" class="data row3 col1">-0.35</td>
      <td id="T__row3_col2" class="data row3 col2">2.94</td>
      <td id="T__row3_col3" class="data row3 col3">0.51</td>
      <td id="T__row3_col4" class="data row3 col4">0.48</td>
      <td id="T__row3_col5" class="data row3 col5">0.83</td>
      <td id="T__row3_col6" class="data row3 col6">0.48</td>
      <td id="T__row3_col7" class="data row3 col7">0.10</td>
      <td id="T__row3_col8" class="data row3 col8">0.90</td>
      <td id="T__row3_col9" class="data row3 col9">0.52</td>
      <td id="T__row3_col10" class="data row3 col10">0.13</td>
      <td id="T__row3_col11" class="data row3 col11">0.21</td>
    </tr>
    <tr>
      <th id="T__level0_row4" class="row_heading level0 row4">CGCNN+P</th>
      <td id="T__row4_col0" class="data row4 col0">0.51</td>
      <td id="T__row4_col1" class="data row4 col1">0.02</td>
      <td id="T__row4_col2" class="data row4 col2">2.38</td>
      <td id="T__row4_col3" class="data row4 col3">0.41</td>
      <td id="T__row4_col4" class="data row4 col4">0.69</td>
      <td id="T__row4_col5" class="data row4 col5">0.78</td>
      <td id="T__row4_col6" class="data row4 col6">0.69</td>
      <td id="T__row4_col7" class="data row4 col7">0.21</td>
      <td id="T__row4_col8" class="data row4 col8">0.79</td>
      <td id="T__row4_col9" class="data row4 col9">0.31</td>
      <td id="T__row4_col10" class="data row4 col10">0.11</td>
      <td id="T__row4_col11" class="data row4 col11">0.18</td>
    </tr>
    <tr>
      <th id="T__level0_row5" class="row_heading level0 row5">CGCNN</th>
      <td id="T__row5_col0" class="data row5 col0">0.52</td>
      <td id="T__row5_col1" class="data row5 col1">-0.61</td>
      <td id="T__row5_col2" class="data row5 col2">2.62</td>
      <td id="T__row5_col3" class="data row5 col3">0.45</td>
      <td id="T__row5_col4" class="data row5 col4">0.60</td>
      <td id="T__row5_col5" class="data row5 col5">0.81</td>
      <td id="T__row5_col6" class="data row5 col6">0.60</td>
      <td id="T__row5_col7" class="data row5 col7">0.15</td>
      <td id="T__row5_col8" class="data row5 col8">0.85</td>
      <td id="T__row5_col9" class="data row5 col9">0.40</td>
      <td id="T__row5_col10" class="data row5 col10">0.14</td>
      <td id="T__row5_col11" class="data row5 col11">0.23</td>
    </tr>
    <tr>
      <th id="T__level0_row6" class="row_heading level0 row6">M3GNet + MEGNet</th>
      <td id="T__row6_col0" class="data row6 col0">0.53</td>
      <td id="T__row6_col1" class="data row6 col1">0.46</td>
      <td id="T__row6_col2" class="data row6 col2">2.65</td>
      <td id="T__row6_col3" class="data row6 col3">0.45</td>
      <td id="T__row6_col4" class="data row6 col4">0.64</td>
      <td id="T__row6_col5" class="data row6 col5">0.80</td>
      <td id="T__row6_col6" class="data row6 col6">0.64</td>
      <td id="T__row6_col7" class="data row6 col7">0.16</td>
      <td id="T__row6_col8" class="data row6 col8">0.84</td>
      <td id="T__row6_col9" class="data row6 col9">0.36</td>
      <td id="T__row6_col10" class="data row6 col10">0.09</td>
      <td id="T__row6_col11" class="data row6 col11">0.13</td>
    </tr>
    <tr>
      <th id="T__level0_row7" class="row_heading level0 row7">M3GNet</th>
      <td id="T__row7_col0" class="data row7 col0">0.58</td>
      <td id="T__row7_col1" class="data row7 col1">0.59</td>
      <td id="T__row7_col2" class="data row7 col2">2.66</td>
      <td id="T__row7_col3" class="data row7 col3">0.45</td>
      <td id="T__row7_col4" class="data row7 col4">0.79</td>
      <td id="T__row7_col5" class="data row7 col5">0.80</td>
      <td id="T__row7_col6" class="data row7 col6">0.79</td>
      <td id="T__row7_col7" class="data row7 col7">0.20</td>
      <td id="T__row7_col8" class="data row7 col8">0.80</td>
      <td id="T__row7_col9" class="data row7 col9">0.21</td>
      <td id="T__row7_col10" class="data row7 col10">0.07</td>
      <td id="T__row7_col11" class="data row7 col11">0.12</td>
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
  #T__row0_col2,
  #T__row0_col3,
  #T__row0_col5,
  #T__row0_col10,
  #T__row1_col7,
  #T__row1_col8,
  #T__row3_col4,
  #T__row3_col6,
  #T__row3_col9,
  #T__row5_col1,
  #T__row5_col11 {
    background-color: #440154;
    color: #f1f1f1;
  }
  #T__row0_col1 {
    background-color: #3c508b;
    color: #f1f1f1;
  }
  #T__row0_col4,
  #T__row0_col6,
  #T__row0_col9 {
    background-color: #472f7d;
    color: #f1f1f1;
  }
  #T__row0_col7,
  #T__row0_col8 {
    background-color: #482374;
    color: #f1f1f1;
  }
  #T__row0_col11 {
    background-color: #424186;
    color: #f1f1f1;
  }
  #T__row1_col0 {
    background-color: #287d8e;
    color: #f1f1f1;
  }
  #T__row1_col1 {
    background-color: #2ab07f;
    color: #f1f1f1;
  }
  #T__row1_col2 {
    background-color: #38598c;
    color: #f1f1f1;
  }
  #T__row1_col3 {
    background-color: #3e4c8a;
    color: #f1f1f1;
  }
  #T__row1_col4,
  #T__row1_col6,
  #T__row1_col9,
  #T__row5_col5,
  #T__row6_col11 {
    background-color: #9bd93c;
    color: #000000;
  }
  #T__row1_col5 {
    background-color: #472a7a;
    color: #f1f1f1;
  }
  #T__row1_col10 {
    background-color: #2a788e;
    color: #f1f1f1;
  }
  #T__row1_col11 {
    background-color: #20a386;
    color: #f1f1f1;
  }
  #T__row2_col0 {
    background-color: #22a785;
    color: #f1f1f1;
  }
  #T__row2_col1 {
    background-color: #238a8d;
    color: #f1f1f1;
  }
  #T__row2_col2,
  #T__row2_col3 {
    background-color: #27808e;
    color: #f1f1f1;
  }
  #T__row2_col4,
  #T__row2_col6,
  #T__row2_col9 {
    background-color: #5cc863;
    color: #000000;
  }
  #T__row2_col5 {
    background-color: #25848e;
    color: #f1f1f1;
  }
  #T__row2_col7,
  #T__row2_col8 {
    background-color: #32648e;
    color: #f1f1f1;
  }
  #T__row2_col10 {
    background-color: #1f978b;
    color: #f1f1f1;
  }
  #T__row2_col11,
  #T__row5_col4,
  #T__row5_col6,
  #T__row5_col9 {
    background-color: #2a778e;
    color: #f1f1f1;
  }
  #T__row3_col0 {
    background-color: #29af7f;
    color: #f1f1f1;
  }
  #T__row3_col1 {
    background-color: #3f4889;
    color: #f1f1f1;
  }
  #T__row3_col2,
  #T__row3_col3,
  #T__row3_col5,
  #T__row3_col7,
  #T__row3_col8,
  #T__row7_col0,
  #T__row7_col1,
  #T__row7_col4,
  #T__row7_col6,
  #T__row7_col9,
  #T__row7_col10,
  #T__row7_col11 {
    background-color: #fde725;
    color: #000000;
  }
  #T__row3_col10 {
    background-color: #46327e;
    color: #f1f1f1;
  }
  #T__row3_col11 {
    background-color: #443983;
    color: #f1f1f1;
  }
  #T__row4_col0 {
    background-color: #52c569;
    color: #000000;
  }
  #T__row4_col1 {
    background-color: #1f968b;
    color: #f1f1f1;
  }
  #T__row4_col2,
  #T__row4_col3 {
    background-color: #25ab82;
    color: #f1f1f1;
  }
  #T__row4_col4,
  #T__row4_col6,
  #T__row4_col9 {
    background-color: #38b977;
    color: #f1f1f1;
  }
  #T__row4_col5 {
    background-color: #3aba76;
    color: #f1f1f1;
  }
  #T__row4_col7,
  #T__row4_col8 {
    background-color: #1f988b;
    color: #f1f1f1;
  }
  #T__row4_col10 {
    background-color: #23898e;
    color: #f1f1f1;
  }
  #T__row4_col11 {
    background-color: #26828e;
    color: #f1f1f1;
  }
  #T__row5_col0 {
    background-color: #56c667;
    color: #000000;
  }
  #T__row5_col2,
  #T__row5_col3,
  #T__row6_col3 {
    background-color: #6ccd5a;
    color: #000000;
  }
  #T__row5_col7,
  #T__row5_col8 {
    background-color: #65cb5e;
    color: #000000;
  }
  #T__row5_col10 {
    background-color: #482475;
    color: #f1f1f1;
  }
  #T__row6_col0,
  #T__row7_col3 {
    background-color: #73d056;
    color: #000000;
  }
  #T__row6_col1 {
    background-color: #b8de29;
    color: #000000;
  }
  #T__row6_col2 {
    background-color: #77d153;
    color: #000000;
  }
  #T__row6_col4,
  #T__row6_col6,
  #T__row6_col9 {
    background-color: #20928c;
    color: #f1f1f1;
  }
  #T__row6_col5 {
    background-color: #98d83e;
    color: #000000;
  }
  #T__row6_col7,
  #T__row6_col8 {
    background-color: #50c46a;
    color: #000000;
  }
  #T__row6_col10 {
    background-color: #7ad151;
    color: #000000;
  }
  #T__row7_col2 {
    background-color: #7fd34e;
    color: #000000;
  }
  #T__row7_col5 {
    background-color: #89d548;
    color: #000000;
  }
  #T__row7_col7,
  #T__row7_col8 {
    background-color: #1fa188;
    color: #f1f1f1;
  }
</style>
