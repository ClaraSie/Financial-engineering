# Financial-engineering 可轉換公司債評價分析
  利用Black & Scholes歐式選擇權、二元樹美式選擇權(tree models)及最小平方蒙地卡羅法(Least-Squares Monte Carlo,LSM)模擬可轉換公司債其演變路徑並進行選擇權評價，找出最適合該可轉換公司債的評價方法。
  研究樣本為聚亨企業股份有限公司國內第四次有擔保轉換公司債，發行日期為民國102年9月30日，發行期間三年，自民國102年9月30日開始發行至民國105年9月30日到期，債券面額為新台幣壹拾萬元，依票面金額發行。

1.讀取發行期間資料為data2022.csv，歷史波動度資料為data2022before.csv，若要更改為自己的資料，請將發行期間與發行前一年分開讀取 </br>
2.在讀取資料時，路徑記得修改為自己電腦的路徑 </br>
3.凍結期間的修改要看T的部分，若開始轉換期間為67，在T的部分則將寫入66，如第8行 </br>
4.最後資料要輸出時，要使用98、99行，99行輸出的路徑也要記得修改</br>
PS.若讀取的資料是excel，記得將副檔名改為csv才能讀取唷
