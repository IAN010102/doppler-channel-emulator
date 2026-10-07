# OFDM 域 ICI 與接收端波束成形：公式與實作對照（`ofdm_ici.js`）

參考：Gopala & Slock, *Doppler Compensation and Beamforming for High Mobility OFDM Transmissions in Multipath*, EURECOM, 2016。
這份文件逐式列出 `ofdm_ici.js` 實作的公式。**編號 O1–O17 是本專案自訂的**；提示詞只提供公式內容，沒有提供文獻的式號，所以「文獻式號」一欄除了「式 24」之外都標為**待確認**。
文獻沒寫明的設定一律是假設（A1、A2 …），見 [ofdm_ici_assumptions.md](ofdm_ici_assumptions.md)。

## 符號與約定（Notation）

| 符號 | 意義 |
|---|---|
| Nr | 接收天線數（均勻線性陣列 ULA，d/λ = 0.5，`a_p(θ) = exp(−j2π p (d/λ) sinθ)`，θ 自法線量起，與模擬器其餘部分相同） |
| N、Δf = fs/N | 子載波數、子載波間隔 |
| L | 路徑數（同一取樣延遲，無延遲擴散 ⇒ 平坦頻率通道） |
| f_i、ε_i = f_i/Δf | 第 i 條路徑的都卜勒偏移（Hz）與正規化都卜勒；`f_i = v cosφ_i fc / c` |
| A_i、θ_i | 路徑複數振幅與到達角 |
| σx²、σn² | 每個子載波的資料功率、雜訊功率（頻域） |
| t | 時間樣本；FFT 窗為 t = 0 … N−1，循環前綴（CP）為 t = −Ncp … −1 |

## 公式

| 編號 | 公式 | 文獻式號 | 程式 |
|---|---|---|---|
| O1 | `Q(x) = (1/N) Σ_{k=0}^{N−1} exp(j2π x k/N) = exp(jπ x (N−1)/N) · sin(πx) / (N sin(πx/N))`（對 x 週期為 N；x ≡ 0 時 Q = 1） | 待確認 | `Qc`、`QbruteForce`（對照用，T23z） |
| O2 | 時域通道 `h(t) = Σ_i A_i a(θ_i) exp(j2π ε_i t/N)`；頻域 `H = Σ_i A_i a(θ_i) Q(ε_i)` | 待確認 | `Hvec`、`makeChannel` |
| O3 | `Y_l = X_l H + Σ_{m≠l} X_m H_ICI(l,m) + W_l`，`H_ICI(l,m) = Σ_i A_i a(θ_i) Q(m − l + ε_i)`（只依 d = (m−l) mod N；d = 0 時就是 H） | 待確認 | `iciVectors` |
| O4 | ICI 協方差 `B_l = Σ_{m≠l, m 被占用} H_ICI(l,m) H_ICI(l,m)ᴴ`（所有子載波皆占用時 B 與 l 無關，省略下標；有保護頻帶時 B_l 隨 l 改變） | 待確認 | `iciCov`、`maskGuard` |
| O5 | `R = σx² B + σn² I`（ICI 加雜訊協方差） | 待確認 | `rMatrix` |
| O6 | 對頻 ε（demodulation）：時間樣本乘 `exp(−j2π ε t/N)`，等價於把 O2–O4 中的 ε_i 換成 ε_i − ε | 待確認 | `iciVectors(ch, eps)`、`demodulate` |
| O7 | 最佳權重 `G_opt = R⁻¹ H` | 待確認 | `optimal` |
| O8 | 最佳 SINR `= σx² Hᴴ R⁻¹ H`；任意權重 `σx² |Gᴴ H|² / (Gᴴ R G)` | 待確認 | `optimal`、`sinrOfWeights` |
| O9 | CP 法：令 `C = Σ_{t=−Ncp}^{−1} Σ_p y_p(t+N) y_p(t)*`，`ε_CP = angle(C)/(2π)`（對頻後 `y'(t+N) y'(t)*` 的相位為零，即結果為實數）。無雜訊的期望值 `E[C] = N σx² Σ_t Σ_p h_p(t+N) h_p(t)*`（CP 使 x(t+N) = x(t)）。適用範圍：`|ε| < 0.5`（相位的 2π 模糊）；多路徑時得到各路徑的相位加權平均；有雜訊時雜訊項期望為零，用同一個相位式估計 | 待確認 | `epsCpAnalytic`、`epsCpEstimate` |
| O10 | 線性近似：`h(t) = h0 + (t − (N−1)/2) h1`（`h0` = 符元中心的通道，`h1` = 該點的斜率；以真實 `h1` 計算，即 genie 版） | 待確認 | `linearModel` |
| O11 | `Ξ_{k,l} = (1/N) Σ_n (n − (N−1)/2) exp(j2π(k−l)n/N)`；ICI 為秩 1（方向 `h1`），`Σ_k |Ξ_{k,l}|² = (N²−1)/12`（Parseval，與 l 無關，`Ξ_{l,l} = 0`） | 待確認 | `linearModel.sumXi2` |
| O12 | `G_lin = (σx² S h1 h1ᴴ + σn² I)⁻¹ h0`，`S = Σ_k |Ξ_{k,l}|²`（提示詞的式子取 σx² = 1） | 待確認 | `linearWeights` |
| O13 | 接收端 FFT：`Y_l = (1/N) Σ_{t=0}^{N−1} y(t) exp(−j2π l t/N)`；傳送 `x(t) = Σ_m X_m exp(j2π m t/N)`；時域雜訊變異數 `N σn²`（使頻域雜訊為 σn²） | 待確認 | `makeSymbol`、`demodulate`、`fft` |
| O14 | `Ĥ = (1/Np) Σ_{l∈導頻} Y_l / X_l`（導頻間距 12，平坦通道對各導頻平均） | 待確認 | `practicalWeights` |
| O15 | `R_yy = (1/Ns) Σ_{l∈S} Y_l Y_lᴴ`（S：用來估協方差的子載波集合；文獻預設 Ns = N；訊號含於 R_yy） | 待確認 | `practicalWeights` |
| O16 | `G_est = R_yy⁻¹ Ĥ`；Ns 可小於 N（連續或均勻間隔選法）；正則化 `none`（Ns ≥ Nr 直接求逆，否則偽逆）或 `dl`（γ = γ_rel σn²） | 待確認 | `practicalWeights` |
| O17 | 比較用接收機：MRC `G = H`（ε = 0，不考慮 ICI）；單天線（ε = 0，只用天線 0）；genie 對頻：`ε* = argmax_ε SINR_opt(ε)`（粗格點加黃金分割），或 `argmax_ε ‖H(ε)‖²`（訊號功率窮舉） | 待確認 | `sinrMrc`（`reproduce_gopala_slock.js` 內同式）、`epsGenie` |
| （式 24） | 簡化的訊號功率最大化近似（選做） | **式 24**；提示詞沒有給出式子，**未實作** | — |

## 驗證（T23）

- T23z：Q 封閉式對暴力加總 < 1e-12；FFT 對 DFT；Parseval `Σ_d |Q(d+ε)|² = 1`；B_l 與 l 無關。
- T23a：N = 64、256，各 4000 個符元，時域逐樣本相乘（與解析公式無關）後做 FFT，量測 ICI 功率與訊號增益，≤ 3 SE。
- T23b：單一路徑、對頻選該路徑的都卜勒 ⇒ ICI 為零，SINR = 陣列增益 Nr × SNR。
- T23c：單一路徑不對頻：ICI 功率 / 總功率 = `1 − |Q(ε)|²`（恆等式），與 T5／UI 的 `1 − sinc²(ε)` 在 N = 1024 的差 < 1e-6（N = 64、256 的差較大，列於輸出）。
- T23d（資訊性）：ICI 協方差最大特徵值的占比。
