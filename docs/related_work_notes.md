# 相關工作筆記（Related work notes）

每篇一行：做了什麼、與本專案差在哪。**閱讀狀態指的是我（Claude）在這個工作階段實際取得的內容**，不是論文本身的完整程度：
- 「讀了全文」＝我讀過論文原文；「只讀摘要」＝我只讀到 arXiv／出版社頁面的摘要（透過網頁抓取工具的轉述，不是原文）；「只看到標題」＝我只有標題或書目資訊；「沒有取得」＝我沒有找到內容。
- 沒讀到的內容我一律不寫。「與本專案差在哪」只在我讀到的範圍內比較；沒有資訊的地方明說。
- 使用者在規格裡標了幾篇的閱讀狀態（Gopala & Slock「讀全文」、Marchese「讀全文」、Feng「只有標題」）；那是使用者的狀態，**我這邊的實際狀態在第一欄另外寫明**，兩者不一定相同。

| 論文 | 我的閱讀狀態 | 做了什麼 | 與本專案差在哪 |
|---|---|---|---|
| Zhang et al., 2011（使用者標示 ChinaCom） | 只看到標題與搜尋結果的一句摘要。搜尋找到的候選是 Ying Zhang et al.，*Multiple Doppler shifts compensation and ICI elimination by beamforming in high-mobility OFDM systems*（2011 年 8 月；IEEE Xplore 文件 6158142）。**這個對應是我的推測，沒有確認它就是使用者指的那篇**；頁面本身我抓不到內容 | （只依搜尋摘要）用陣列波束成形與干擾抑制補償高速 OFDM 中的多重都卜勒偏移、消除 ICI | 本專案舊的 ICI 模型把 ICI 當功率項、接收端不做對頻，沒有涵蓋這類處理；`ofdm_ici.js` 的最佳合併與 AD-MF 是同類想法的簡化實作，但我沒有讀過這篇，無法比較方法細節 |
| Gopala & Slock, 2016（EURECOM；*Doppler Compensation and Beamforming for High Mobility OFDM Transmissions in Multipath*） | **使用者標示讀了全文；我沒有讀過原文**，只看到提示詞摘錄的公式與模擬設定。另外搜尋找到同名的 Springer 章節（`10.1007/978-3-319-40352-6_24`），它與 EURECOM 版本是否相同我沒有檢查 | （依提示詞摘錄）SIMO OFDM、多路徑各有自己的都卜勒；頻域模型含 ICI；接收端最佳波束成形 `G = R⁻¹H`；對頻選法（不對頻、LoS、CP 法、genie）；導頻估計 Ĥ 與樣本 R_yy 的實務接收機；通道線性變化的近似版 | 本專案的 `ofdm_ici.js` 實作這些公式（編號 O1–O17，文獻式號待確認），並在其上加了：獨立的時域驗證（T23）、假設與敏感度（`docs/ofdm_ici_assumptions.md`）、Estimated BF 損失的分解（含訊號的協方差造成自我抵消）、導頻殘差協方差估計器（O21）。原文的圖我沒看到，所以沒有逐圖比對 |
| Guo et al., 2017, IEEE TVT（*High-Mobility OFDM Downlink Transmission With Large-Scale Antenna Array*，Wei Guo、Weile Zhang、Pengcheng Mu、Feifei Gao；arXiv 1701.03221） | 只讀摘要（arXiv 摘要頁） | 針對大規模天線陣列的高移動 OFDM 下行：用波束成形網路在空間域分離多個載波頻偏（振盪器頻偏加多個都卜勒偏移），聯合估計都卜勒與振盪器頻偏，再對每個波束分支做傳統的單一頻偏補償與通道估計；不需要複雜的時變通道估計 | 本專案的 AD-MF 是這個「每條路徑一個分支」想法的 **genie 版**（真實角度與都卜勒、不估計、不做頻偏補償以外的處理，分支用 MRC 合併）；該文有估計器（依摘要），本專案沒有。除此之外細節我沒讀 |
| Ge et al.（arXiv 1809.00137，*Beamforming Network Optimization for Reducing Channel Time Variation in High-Mobility Massive MIMO*；Yinghao Ge、Weile Zhang、Feifei Gao、Shun Zhang、Xiaoli Ma） | 只讀摘要 | 高鐵上行到基地台的快速時變通道：推導通道功率譜密度，用大型均勻線性陣列的（發射）波束成形分離都卜勒，並引入一個共同可調的振幅與相位參數來最佳化波束成形網路（放寬傳統匹配濾波波束成形器的限制），給出使都卜勒擴展最小的最佳參數閉式解 | 它設計的是降低通道時變的波束成形網路（摘要中是發射端波束成形、上行），本專案做的是接收端合併與協方差估計；本專案沒有實作它的參數最佳化。其餘我沒讀 |
| Yan et al., 2024, IET Communications | 只看到標題與搜尋摘要（出版社頁面回應 403，抓不到）。搜尋找到的候選是 *Doppler spread analysis for suppressing channel time variation in high-mobility massive MIMO V2V communications*（作者名在搜尋摘要中是 Zeyu Yan）；**這個對應是我的推測，沒有確認** | （只依搜尋摘要）高移動大規模 MIMO 的 V2V 通訊：接收端多個都卜勒頻偏的疊加使衰落通道出現時變（時間選擇性衰落），分析都卜勒擴展以抑制通道時變 | 本專案的通道是平坦頻率、直線軌道的單一實現模型，沒有 V2V 的雙端移動；也沒有做都卜勒擴展的解析分析。其餘我沒讀 |
| Feng et al., 2022, IEEE TCOM | **沒有取得**：使用者標示「只有標題」，但標題也沒有提供給我，我搜尋也沒有找到這篇，所以連標題都沒有 | 無法描述 | 無法比較；待使用者提供標題或全文後再補 |
| Marchese et al., 2026（arXiv 2601.12970，*6G OFDM Communications with High Mobility Transceivers and Scatterers via Angle-Domain Processing and Deep Learning*；Mauro Marchese、Musa Furkan Keskin、Henk Wymeersch、Pietro Savazzi；2026 年 1 月 19 日提交） | **使用者標示讀了全文；我只讀到摘要**（arXiv 摘要頁） | OFDM 接收機用角度域處理在高移動情境分離多徑，結合導頻估計 DoA、延遲、通道增益，再用決策導向的細化做都卜勒估計；比較兩種初始都卜勒估計（基於 EVM 與深度學習），模擬指出深度學習版在到 1000 km/h 仍維持穩定的效能（依摘要） | 本專案的 AD-MF 是同一類角度域處理的最簡 genie 基準（真實角度與都卜勒，沒有 DoA、延遲、增益或都卜勒的估計，沒有決策導向與深度學習，沒有延遲擴散）。本專案沒有重現該文的任何結果 |
| *Minimum Eigenvalue Based Covariance Matrix Estimation with Limited Samples*（arXiv 2306.11332；Jing Qian、Juening Jin、Hao Wang，2023 年 6 月 20 日提交） | 只讀摘要 | 針對樣本數有限時干擾抑制合併（IRC）接收機的協方差估計：不用最小化均方誤差，改以最小特徵值準則最佳化正則化參數，直接針對 IRC 接收機效能，並推導較易算的最小特徵值下界；模擬顯示在互資訊指標上可接近 oracle 估計器（依摘要） | 本專案的導頻殘差估計器用固定的對角加載（γ_rel）或 Ledoit–Wolf 收縮，沒有做最小特徵值準則的正則化最佳化；兩者的比較沒有做。其餘我沒讀 |

## 備註

- 上表「做了什麼」中標「依搜尋摘要」或「依摘要」的內容，都只是摘要層級；要引用請先讀原文。
- 本專案在這一輪另外引用的 Ledoit & Wolf (2004)（J. Multivariate Analysis 88(2): 365–411）只用於收縮公式（`ofdm_ici.js` 的 `ledoitWolfShrink`），這份筆記沒有把它列入，因為不是 OFDM 的相關工作。
- 網址（arXiv 摘要頁）：<https://arxiv.org/abs/1809.00137>、<https://arxiv.org/abs/2601.12970>、<https://arxiv.org/abs/2306.11332>、<https://arxiv.org/abs/1701.03221>。
