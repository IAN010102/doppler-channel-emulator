# 最終報告：陣列訊號處理 × 高移動性都卜勒模擬器（分支 `unified-doppler`）

Written for: 專案負責人（作為交付與論文引用的依據）。狀態以 HEAD（`7fff66a` 之後的最終 commit）為準；未 merge、未部署、未推送。
**MATLAB 部分未在本環境執行**（無 MATLAB）；網頁端對 MATLAB 資料的自洽性只在 Node 端驗證（T20）。

---

## 1. 各 commit 的變更與關鍵函式

| commit | 內容與關鍵函式 |
|---|---|
| 14 `e4a3e49` | MMSE 分成 MMSE-M（模型式，r_xd = P_s a）與 MMSE-P（導頻式，`rxdHat`），選單／掃描／CSV；T11 |
| 15 `a399d1f` | `smiSingular`（pinv 預設｜clamp 舊）、`pinvHermitian`、秩與條件數讀數、奇異提示；T12 |
| 16 `941a301` | 目標符元隨調變、`Core.symbolLevel` / `Sys.symbolEvm`、符元層 EVM／SER；T13 |
| 17 `8f959a6` | 預設改為 unified + signalFree（舊模型標為 v4 analytic）、golden 快照 T14 |
| 18 `d3a4792` | 只做診斷（`tests/diag_c18_mmsep.js`），**未實作** MMSE-P 頻偏補償 |
| 19 `6692edd` | `wienerSolve`（L<N 用 pinv）、相對閾值 `invertMatrix`（舊規則 `invertMatrixAbs`）、`jamWave`、CSV `jam_wave`、golden 重產；T16、T13c 閘門 |
| 20 `1c22615` | `pointErrDeg`（`steer(θ̂₁+δθ)`；MMSE-P 不受影響）、讀數、CSV；T17 |
| 21 `ad0f369` | DOA 空間譜 `doaSpectra`（Capon＋pinv、MUSIC 2 源）、`spectrumPeaks`、`Sys.computeDoa`；T18 |
| 22 `e0a09c1` | `experiments.js`（E1–E5）、`i18n.js`（繁中／英文、提示 ≤40 字）、`d_min` 控制；T19 |
| 22b `4eacab5` | `covSource`（`theoryCov`）、原始 EVM `evmRaw`、`angleSource` music、實驗 E0、角度慣例文件；T21；Lang 觀察器無窮迴圈修正 |
| 23 `0567c16` | `data/matlab_cases.json`（23 案例）、`curve_*.csv`、`README_data.md` |
| 24 `2b741bd` | MATLAB：`verify_cases`、`mc_independent`、`plot_figures`、`lab_mmse_mvdr_demo`、五關學習階梯＋`check_my_work`＋`selftest`；T20 |
| 25 `0ca1788` | 真瀏覽器驗證（Edge＋playwright-core）與截圖 |
| 22c `c85007b` | `dEffT/dEffJ`（δθ_eff）、`pointingMode`、實驗 E6、CSV `pointing_mode`、`delta_theta_eff_deg`；T22 |
| `384b670` | `.kv .k` 換行修正、瀏覽器驗證重跑、診斷 B1–B4（`docs/diagnostics/`）、T17b2（資訊性）、golden CHANGELOG 標題 |
| `64ed59e`、`24e0e60` | 敏感度診斷 S1–S3、E6 文字註明 (fc, Δf, d_min, τ) 與適用範圍、PARAMS.md §21 定義 |
| `7fff66a` | T11c／T13a／T17b 判準修訂（見第 2 項） |

## 2. 完整測試結果（`node tests/run_all.js`，最終狀態）

**198 / 198 項 PASS，0 FAIL。**（資訊性項目不計入。）

| 測試 | PASS | FAIL | 時間 | 測試 | PASS | FAIL | 時間 |
|---|---|---|---|---|---|---|---|
| T0a 與舊頁面逐位相同 | 1 | 0 | 1 s | T12 小快照 pinv | 15 | 0 | 257 s |
| T0b UI 煙霧測試 | 1 | 0 | 7 s | T13 符元層 EVM | 11 | 0 | 100 s |
| T1 unified vs legacy | 6 | 0 | 207 s | T14 golden | 13 | 0 | 0.2 s |
| T2 單徑不變性 | 5 | 0 | 42 s | T16 奇異／相對閾值 | 8 | 0 | 74 s |
| T3 老化 | 10 | 0 | 816 s | T17 指向偏差 | 10 | 0 | 1261 s |
| T4 擴散 ICI | 12 | 0 | 0.1 s | T18 DOA | 5 | 0 | 15 s |
| T5 時域 OFDM | 7 | 0 | 1 s | T19 文字／實驗 | 36 | 0 | 376 s |
| T6 DL γ_rel（僅資訊） | — | — | 171 s | T20 MATLAB 資料 | 17 | 0 | 2 s |
| T7 幾何 | 4 | 0 | 0 s | T21 講義基準 | 6 | 0 | 154 s |
| T8 積分相位 | 5 | 0 | 0.5 s | T22 移動性指向 | 6 | 0 | 38 s |
| T9 numerology | 7 | 0 | 0.4 s | T10 種子／CSV | 8 | 0 | 169 s |
| T11 MMSE | 5 | 0 | 47 s | | | | |

**這輪沒有未通過項目。但要誠實說明：先前的 4 項失敗（T11c、T13a、T17b×2）不是因為程式行為改變而通過，而是因為判準依調查結果修訂。** 倍數（SE 的 1 倍、3 倍）都沒有放寬：

| 項目 | 修改前 | 修改後 | 依據 |
|---|---|---|---|
| T11c | L = 4, 8, 12, 24, 48, 100 全部檢定單調上升（4→8 下降 7.7 dB 而失敗） | 只檢定 L ≥ N+2（12, 24, 48, 100），並輸出 L = 7, 8, 9 的 SINR、κ、‖w‖ | B1：L=N 是低谷（κ 中位 5.5e5 對 L=9 的 1.9e5；‖w‖ 1.12 對 0.75）。**與你的指示有一處差異**：只排除 \|L−N\|≤1 仍會留下 L=4，而 L=4 本身就在低谷左側（B1：L=4→7 由 25.9 降到 21.9 dB），4→12 仍會下降，所以 L=4 也一併排除 |
| T13a | 單一區塊 300 次，\|mean\| ≤ 3 SE | 10 個獨立區塊（各 300 次）以 Stouffer 合併 Z，\|Z\| ≤ 3；輸出各區塊 z 分佈，sd > 1.5 則註記 | signalFree v=300 區塊 0 的 z=3.72，其餘 9 塊 \|z\| ≤ 1.45（B2）。合併後 Z=1.31；各格區塊 z 的 sd 為 0.70 / 1.40 / 0.70 / 0.60，無一超過 1.5，**沒有 SE 偏低的證據**（SE 為樣本標準差／√300，無重尾修正、無 bootstrap） |
| T17b（signalFree，δθ=0 與 1°） | 與母體解比較 | 與「母體解 + Reed–Mallett–Brennan 期望損失 −0.00152 dB」比較 | T17b2：與理論差 −0.0007 dB（z = −2.07）、−0.0011 dB（z = −1.78），在 3 SE 內；對 0 比較原本是 7.3 與 4.3 SE。δθ ≥ 3° 與全部 withSignal 維持原判準與原結果 |

仍屬觀察而未動的項目（資訊性）：T16c 的 L=8→9 步階（MMSE-P +4.08 ± 0.16 dB；MMSE-M +0.57；SMI signalFree +1.74），原因見 `docs/diagnostics/b1_mmsep_L.txt`（L=N 低谷）。

**資訊性項目**（不影響 PASS/FAIL）：T1b、T3c、T6、T11c 觀察、T12 的 SINR–L、T13c 低錯誤數列、T16b/c、T17b2/d/e、T18c、T21d、T22c，輸出在 `node tests/run_all.js` 的 INFORMATIONAL 區。

## 3. Commit 18 診斷表（B4；withSignal 列為修正版）

設定：unified、L = 100、SNR 20、SIR −10、θ₂ = 40°、QPSK、各格 500 次，SINR 單位 dB ± SE。MMSE-P 的 R̂ 一律含目標，因此 signalFree 與 withSignal 完全相同，下表只列一份。
**更正：**舊診斷腳本在 withSignal 時把目標算了兩次（`rr` 已含目標又加 `tr`），當時 withSignal 列的 genie 與 lag-1 數字作廢；已修正（`diag_c18_mmsep.js`、`diag_b4_afc.js`）。

欄位：(i) 現有 MMSE-P、(ii) genie 去旋轉（用真實 LoS 積分相位）、cycles = 視窗內 LoS 旋轉圈數、ω = 每快照平均相位增量（rad，皆 < π）、(iv) 規格的 lag-1 估計（ω̂ 誤差為 RMS，rad）、(v) 週期圖法（網格 π/(8L)＋拋物線內插）。

```
K = 20 dB, theta1 = 0, trainMode = signalFree   (L = 100, 500 trials; SINR dB ± SE; omega in rad/snapshot)
v km/h   cycles    omega    (i) MMSE-P  gap opt    (ii) genie  gap opt    (iv) lag-1  gap opt   RMS err   (v) period.  gap opt   RMS err
0          0.00   0.0000  28.61 ± 0.02    -0.32  28.61 ± 0.02    -0.32 -15.30 ± 0.48   -44.23   8.37e-1  26.68 ± 0.31    -2.25   2.50e-1
30         0.98   0.0623 -20.02 ± 0.28   -48.95  28.61 ± 0.02    -0.32 -14.86 ± 0.45   -43.79   8.38e-1  26.66 ± 0.29    -2.27   2.84e-1
100        3.27   0.2075  -9.65 ± 0.21   -38.58  28.59 ± 0.02    -0.34 -12.01 ± 0.40   -40.94   8.34e-1  26.27 ± 0.29    -2.67   2.55e-1
300        9.81   0.6225 -12.36 ± 0.26   -41.29  28.51 ± 0.02    -0.42  -8.72 ± 0.36   -37.65   8.43e-1  25.63 ± 0.27    -3.30   2.45e-1
```

```
K = 20 dB, theta1 = 45, trainMode = signalFree   (L = 100, 500 trials; SINR dB ± SE; omega in rad/snapshot)
v km/h   cycles    omega    (i) MMSE-P  gap opt    (ii) genie  gap opt    (iv) lag-1  gap opt   RMS err   (v) period.  gap opt   RMS err
0          0.00   0.0000  21.61 ± 0.04    -0.32  21.61 ± 0.04    -0.32 -15.32 ± 0.47   -37.26   9.24e-1  15.58 ± 0.57    -6.35   5.61e-1
30         0.69   0.0440   1.36 ± 0.09   -20.62  21.38 ± 0.03    -0.60 -12.21 ± 0.39   -34.19   9.27e-1  14.72 ± 0.54    -7.26   5.91e-1
100        2.32   0.1470  -3.83 ± 0.18   -25.78  19.45 ± 0.04    -2.50  -9.61 ± 0.39   -31.56   9.35e-1  13.61 ± 0.49    -8.34   5.97e-1
300        6.97   0.4423 -13.97 ± 0.28   -35.94  14.76 ± 0.05    -7.22 -10.20 ± 0.35   -32.17   9.47e-1  10.71 ± 0.41   -11.27   6.14e-1
```

```
K = inf (400 dB), theta1 = 0, trainMode = signalFree   (L = 100, 500 trials; SINR dB ± SE; omega in rad/snapshot)
v km/h   cycles    omega    (i) MMSE-P  gap opt    (ii) genie  gap opt    (iv) lag-1  gap opt   RMS err   (v) period.  gap opt   RMS err
0          0.00   0.0000  28.62 ± 0.01    -0.32  28.62 ± 0.01    -0.32 -15.22 ± 0.47   -44.16   8.36e-1  26.96 ± 0.26    -1.98   2.19e-1
30         0.98   0.0623 -21.83 ± 0.24   -50.77  28.62 ± 0.01    -0.32 -15.70 ± 0.46   -44.64   8.37e-1  27.08 ± 0.24    -1.86   2.23e-1
100        3.27   0.2075 -10.88 ± 0.13   -39.82  28.60 ± 0.01    -0.34 -15.44 ± 0.47   -44.38   8.33e-1  27.10 ± 0.26    -1.84   2.34e-1
300        9.81   0.6225 -22.31 ± 0.24   -51.25  28.56 ± 0.01    -0.38 -16.10 ± 0.46   -45.04   8.42e-1  26.93 ± 0.21    -2.01   2.15e-1
```

```
K = inf (400 dB), theta1 = 45, trainMode = signalFree   (L = 100, 500 trials; SINR dB ± SE; omega in rad/snapshot)
v km/h   cycles    omega    (i) MMSE-P  gap opt    (ii) genie  gap opt    (iv) lag-1  gap opt   RMS err   (v) period.  gap opt   RMS err
0          0.00   0.0000  21.65 ± 0.01    -0.32  21.65 ± 0.01    -0.32 -15.08 ± 0.47   -37.06   9.34e-1  15.70 ± 0.56    -6.28   5.79e-1
30         0.69   0.0440   3.35 ± 0.08   -18.63  21.51 ± 0.01    -0.46 -14.99 ± 0.45   -36.96   9.39e-1  15.42 ± 0.57    -6.56   5.73e-1
100        2.32   0.1470  -6.89 ± 0.10   -28.86  20.22 ± 0.01    -1.76 -14.66 ± 0.45   -36.64   9.44e-1  14.26 ± 0.56    -7.71   5.84e-1
300        6.97   0.4423 -20.00 ± 0.23   -41.98  15.48 ± 0.02    -6.50 -14.87 ± 0.41   -36.84   9.50e-1  11.17 ± 0.47   -10.81   5.99e-1
```

```
Criteria of T15 evaluated with the periodogram estimator (K = 20 dB, L = 100, signalFree = withSignal for MMSE-P)
theta1 = 0: T15a  v = 0, on - off = -1.93 dB (SE 0.31, |z| = 6.2) -> would FAIL;  T15b  v = 300: on - off = 37.99 dB (>= 10: yes), genie - on = 2.87 dB (< 3: yes)
theta1 = 45: T15a  v = 0, on - off = -6.03 dB (SE 0.56, |z| = 10.7) -> would FAIL;  T15b  v = 300: on - off = 24.67 dB (>= 10: yes), genie - on = 4.05 dB (< 3: no)
```

**結論（不實作 AFC）：**
- 假設「高速崩潰是因為 r̂_xd 平均時沒有補償公共相位旋轉」**被支持**：genie 去旋轉後 θ₁=0° 幾乎回到 SINR_opt（差 0.3–0.4 dB）；θ₁=45° 在 300 km/h 仍差 6.5–7.2 dB（剩餘來自 LoS 以外的機制，如角度漂移與擴散）。
- 規格的 lag-1 估計器不可用（ω̂ RMS 誤差 0.83–0.95 rad，SIR −10 dB 下幾乎隨機，SINR −8 到 −16 dB）。
- 週期圖法在 v ≥ 30 km/h 大幅改善（θ₁=0°、K=20 dB、300 km/h：25.6 dB，比不補償高 38 dB，距 genie 2.87 dB），但依先前 T15 判準：**T15a（v=0，on 與 off 差在 3 SE 內）不過**（θ₁=0°：−1.93 dB、|z|=6.2；θ₁=45°：−6.03 dB、|z|=10.7）；**T15b 在 θ₁=45° 不過**（與 genie 差 4.05 dB，> 3 dB）。原因：SIR −10 dB 下約 1.0–1.3 % 的試驗峰值落在錯誤頻率（v=0 時誤差中位數 1.6e-3 rad、第 90 百分位 5.5e-3 rad，但有離群值）。因此停在診斷，**MMSE-P 在本專案中標示為「無頻偏補償」**。

## 4. golden 變動清單與原因

- Commit 17 建立快照（預設切為 unified + signalFree）；Commit 19 是唯一一次改動數值的 commit。
- **Commit 19**：`jamWave` 預設改為 `'gaussian'`，干擾符元抽取不同亂數；FOURIER 權重與干擾波形無關，SINR 不變。逐種子 SINR（dB，種子 11/22/33/44，舊 → 新；EVM、SINR_opt、SER 隨之改變）：

```
default(signalFree)/FOURIER  13.23->13.23  15.75->15.75  17.70->17.70  13.38->13.38
default(signalFree)/MMSE     -2.64->-3.50  -4.95->-4.29  0.00->-0.83  -9.43->-7.49
default(signalFree)/MMSEP    -8.33->-11.20  -0.68->-7.61  -3.44->-6.96  -5.51->-10.81
default(signalFree)/SMI      27.81->28.35  28.23->28.19  28.57->28.42  28.76->29.13
default(signalFree)/DL       28.30->28.44  28.57->28.24  28.92->28.70  28.95->29.26
default(signalFree)/BEAMSPACE 27.11->27.14  27.47->27.04  27.33->27.15  26.94->27.34
withSignal/FOURIER           13.23->13.23  15.75->15.75  17.70->17.70  13.38->13.38
withSignal/MMSE              -2.64->-3.50  -4.95->-4.29  0.00->-0.83  -9.43->-7.49
withSignal/MMSEP             -8.33->-11.20  -0.68->-7.61  -3.44->-6.96  -5.51->-10.81
withSignal/SMI               -2.64->-3.50  -4.95->-4.29  0.00->-0.83  -9.43->-7.49
withSignal/DL                13.55->13.56  9.64->10.98  14.32->15.12  8.85->8.72
withSignal/BEAMSPACE         5.27->6.17  1.19->2.42  6.82->5.28  -2.54->-3.14
```
- Commit 20–22c 與其後所有 commit：golden **未變**（T14、T17a 逐位比對通過）。未實作 AFC，因此沒有 AFC 造成的 golden 變動。

## 5. 資訊性結果摘要

- **T17d（指向偏差，L=100）**：signalFree 時 SMI 28.6→28.0（3°）→26.8（5°）、DL 28.9→28.2→27.1、BEAMSPACE 28.8→28.7→28.5，MMSE-P 不受影響（28.6）。withSignal 時 SMI −5.5→−19.7→−21.5、DL 13.9→−0.3→−5.4、BEAMSPACE 6.9→−13.0→−17.4。K=40 dB 時 signalFree 與 withSignal 的差在 δθ=0° 就已超過 3 dB。
- **T17e**：withSignal、δθ=3° 的 DL，最佳 γ_rel = +28 dB（20.51 dB），預設 +10 dB 為 −0.28 dB（預設值未改）。
- **T18c（DOA 峰值誤差）**：L ≥ 12、分離 ≥ 20° 時 MUSIC 與 Capon 目標角誤差都約 0.2°（網格 0.5°，量化極限 0.25°）；分離 5° 時 MUSIC 1.29°、Capon 1.85°；L=4 時 Capon 0.37°、MUSIC 0.23°。
- **T21d（angleSource=music，SINR 名義角／MUSIC 角）**：signalFree 的 SMI 從 L=4 的 22.4 升到 L=100 的 28.6 dB（L=8 低谷 21.7）；signalFree 的 SMI／DL／BEAMSPACE 兩種角度相差 ≤ 0.1 dB；含訊號訓練時 MUSIC 角較高（SMI／MMSE-M 最多約 +1.0 dB，DL 約 +0.1 到 +0.7 dB）；MUSIC 角誤差中位數 0.00°、第 90 百分位 0.50°。
- **T6（DL γ_rel）**：signalFree 時最佳 γ_rel 約 +10 到 +12 dB；withSignal 時 L=100 最佳 +26 dB（21.18 dB），L=12 到範圍上限都在上升。
- **E6 與敏感度**：見第 6 項；S1（240 格）沒有任何一格老化損失 ≥ ICI 損失；S2：300 km/h 時 DL 損失 4.91 dB、SMI 1.35 dB，DL 零陷淺約 9 dB、寬度相近；S3：withSignal 的上升在 K→∞ 消失。
- **B1–B3**：見 `docs/diagnostics/README.md`。

## 6. 各實驗預設（E0–E6）的實測數字與最終文字

所有數字由 `tests/experiment_measure.js`（固定種子）量測，存於 `data/experiment_numbers.json`，T19d 重新量測並比對。**實測數字（`data/experiment_numbers.json`）**

```
{
 "E0": {
  "MMSE-M (theory R)": {
   "sinr": 38.97,
   "opt": 38.97,
   "evmRaw_pct": 1.13,
   "bias": 0.000127
  },
  "MMSE-M (sample R)": {
   "sinr": 21.88,
   "opt": 38.97,
   "evmRaw_pct": 8.29,
   "bias": 0.006781
  },
  "SMI withSignal (theory R) = MVDR B": {
   "sinr": 38.97,
   "opt": 38.97,
   "evmRaw_pct": 1.13,
   "bias": 0
  },
  "SMI withSignal (sample R)": {
   "sinr": 21.88,
   "opt": 38.97,
   "evmRaw_pct": 8.2,
   "bias": 0
  },
  "SMI signalFree (theory R) = MVDR A": {
   "sinr": 38.97,
   "opt": 38.97,
   "evmRaw_pct": 1.13,
   "bias": 0
  },
  "SMI signalFree (sample R)": {
   "sinr": 38.94,
   "opt": 38.97,
   "evmRaw_pct": 1.13,
   "bias": 0
  }
 },
 "E1_28": {
  "byV": {
   "0": {
    "evm_pct": 3.75,
    "ici_pct": 0
   },
   "100": {
    "evm_pct": 5.43,
    "ici_pct": 3.92
   },
   "200": {
    "evm_pct": 8.68,
    "ici_pct": 7.83
   },
   "300": {
    "evm_pct": 12.31,
    "ici_pct": 11.72
   }
  },
  "crossing_speed_kmh": null
 },
 "E1_5": {
  "byV": {
   "0": {
    "evm_pct": 3.75,
    "ici_pct": 0
   },
   "100": {
    "evm_pct": 6.75,
    "ici_pct": 5.59
   },
   "200": {
    "evm_pct": 11.79,
    "ici_pct": 11.17
   },
   "300": {
    "evm_pct": 17.12,
    "ici_pct": 16.7
   }
  },
  "crossing_speed_kmh": 230
 },
 "E2": {
  "SMI": {
   "sinr": 22.46,
   "rank": 4
  },
  "DL": {
   "sinr": 27.81,
   "rank": 4
  },
  "BEAMSPACE": {
   "sinr": 26.36,
   "rank": 4
  }
 },
 "E3": {
  "ref_signalFree_SMI_d0": 28.57,
  "SMI": {
   "d0": -5.56,
   "d3": -19.75
  },
  "DL": {
   "d0": 13.98,
   "d3": -0.27
  },
  "MMSEP": {
   "d0": 28.58,
   "d3": 28.58
  }
 },
 "E4": {
  "signalFree": {
   "0": 28.59,
   "1": 28.52,
   "3": 27.96,
   "5": 26.8
  },
  "withSignal": {
   "0": -5.72,
   "1": -12.03,
   "3": -19.89,
   "5": -21.61
  }
 },
 "E5": {
  "d30": {
   "SMI": {
    "0": 28.31,
    "5": 27.07,
    "10": 25.48
   },
   "DL": {
    "0": 28.31,
    "5": 26.11,
    "10": 23.75
   }
  },
  "d5": {
   "SMI": {
    "0": 28.13,
    "5": 27.37,
    "10": 25.87
   },
   "DL": {
    "0": 26.19,
    "5": 21,
    "10": 17.17
   }
  }
 },
 "E6": {
  "byV": {
   "0": {
    "dEffT": 0,
    "dEffJ": 0,
    "ratio_pct": 0,
    "bw3": 13.84,
    "iciDb": -300,
    "agingDb": -300,
    "smiSF": 28.61,
    "smiWS": -5.61,
    "dlSF": 28.88,
    "mmseP": 28.63
   },
   "50": {
    "dEffT": 0.09,
    "dEffJ": -0.2,
    "ratio_pct": 0.67,
    "bw3": 13.83,
    "iciDb": -31.6,
    "agingDb": -39.21,
    "smiSF": 28.24,
    "smiWS": 4.62,
    "dlSF": 28.62,
    "mmseP": -3.47
   },
   "100": {
    "dEffT": 0.18,
    "dEffJ": -0.39,
    "ratio_pct": 1.33,
    "bw3": 13.86,
    "iciDb": -25.58,
    "agingDb": -35.08,
    "smiSF": 27.73,
    "smiWS": 6.99,
    "dlSF": 27.81,
    "mmseP": -5.21
   },
   "200": {
    "dEffT": 0.37,
    "dEffJ": -0.78,
    "ratio_pct": 2.63,
    "bw3": 13.92,
    "iciDb": -19.58,
    "agingDb": -32.94,
    "smiSF": 27.24,
    "smiWS": 5.9,
    "dlSF": 25.61,
    "mmseP": -8.35
   },
   "300": {
    "dEffT": 0.54,
    "dEffJ": -1.15,
    "ratio_pct": 3.89,
    "bw3": 14,
    "iciDb": -16.08,
    "agingDb": -33.38,
    "smiSF": 27.36,
    "smiWS": 5.05,
    "dlSF": 23.91,
    "mmseP": -11.17
   }
  },
  "cross_ici_over_aging_kmh": 25,
  "peak_smiWS": {
   "v": 175,
   "sinr": 7.14
  },
  "tau0": {
   "smiSF_300": 28.51,
   "smiWS_300": 6.58,
   "smiWS_100": 0.11
  },
  "scope": {
   "cells": 240,
   "ge": 0,
   "maxRatio": 0.56,
   "at": {
    "fc": 28,
    "scs": 120,
    "dmin": 10,
    "tau": 10,
    "v": 300,
    "algo": "SMI",
    "aging": 5.53,
    "ici": 9.91
   }
  }
 }
}
```

E6 使用 fc = 5 GHz、Δf = 15 kHz、d_min = 10 m、τ = 10 ms、θ₁ = 20°、θ₂ = −30°、K = 20 dB、L = 100、N = 8。


#### E0 [zh] E0 講義基準
- look: 用講義的理論協方差時，MMSE-M、MVDR 方法 B（含目標的 R_r）、方法 A（不含目標的 R_u）的 SINR 都是 38.97 dB，等於 SINR_opt（38.97 dB），而且 MMSE 與方法 B 的權重方向相同。改用 1000 筆快照的樣本協方差後，含目標的 MMSE-M 與 SMI 只剩 21.88 dB，不含目標的 SMI 仍有 38.94 dB。原始 EVM（不做增益正規化）約 1.13%；MMSE 的輸出縮放偏差 |g′−1| 平均只有 0.000127，MVDR 為 0。
- why: 理論協方差沒有抽樣誤差，MMSE 的 R_r⁻¹a 與 MVDR-B 的 R_r⁻¹a/(aᴴR_r⁻¹a) 只差一個純量，所以方向相同，也都等於最佳權重。用樣本估計含目標的協方差時，估計誤差會被當成可以壓掉的干擾（自我抵消）；不含目標的訓練資料沒有這個問題。
- real: 講義的 MMSE 程式用理論值，所以結果漂亮；真實接收機只有樣本協方差，若訓練資料含訊號就會自我抵消。可以先拿這個基準對答案（MATLAB 的 lab_mmse_mvdr_demo.m），再比較樣本版本差多少。

#### E0 [en] E0 Lecture baseline
- look: With the theoretical covariance of the lecture, MMSE-M, MVDR method B (R_r with the target) and method A (R_u without it) all reach 38.97 dB, equal to SINR_opt (38.97 dB), and MMSE and method B have the same weight direction. With the sample covariance of 1000 snapshots, MMSE-M and SMI that train with the target drop to 21.88 dB, while SMI without the target keeps 38.94 dB. The raw EVM (no gain normalisation) is about 1.13%; the output scaling bias |g′−1| of MMSE is only 0.000127 on average, that of MVDR is 0.
- why: The theoretical covariance has no sampling error; the MMSE weight R_r⁻¹a and the MVDR-B weight R_r⁻¹a/(aᴴR_r⁻¹a) differ only by a scalar, so they point the same way and both equal the optimum. When the covariance that contains the target is estimated from samples, the estimation error is treated as interference that can be suppressed (self-nulling); training data without the target do not have this problem.
- real: The lecture MMSE program uses theoretical values, hence the clean result; a real receiver has only the sample covariance, and with the signal in the training data it suffers self-nulling. Check your program against this baseline first (lab_mmse_mvdr_demo.m in MATLAB), then see how far the sample version falls short.

#### E1 [zh] E1 高速鐵路
- look: 把車速從 0 拉到 300 km/h，EVM 曲線一路上升並貼近虛線（ICI 極限）。300 km/h 時：28 GHz / 120 kHz 的 EVM 約 12.31%（ICI 極限 11.72%），5 GHz / 15 kHz 約 17.12%（極限 16.7%）。
- why: 車速越快，都卜勒頻移越大，子載波不再正交，產生載波間干擾（ICI）。天線陣列只能濾掉來自其他方向的干擾，濾不掉這種自己造成的干擾，所以曲線會貼住 ICI 極限。
- real: ICI 由頻移與子載波間隔的比值決定：28 GHz / 120 kHz 的 fc/Δf 比 5 GHz / 15 kHz 小，所以 300 km/h 仍低於 16-QAM 的 13.1% 門檻，而 5 GHz / 15 kHz 在約 230 km/h 就超過了。

#### E1 [en] E1 High-speed rail
- look: Raise the speed from 0 to 300 km/h: the EVM curve rises and hugs the dashed line (the ICI limit). At 300 km/h the EVM is about 12.31% at 28 GHz / 120 kHz (ICI limit 11.72%) and about 17.12% at 5 GHz / 15 kHz (limit 16.7%).
- why: The faster the vehicle, the larger the Doppler shift; the subcarriers are no longer orthogonal and inter-carrier interference (ICI) appears. An antenna array can only remove interference from other directions, not this self-made one, so the curve sticks to the ICI limit.
- real: ICI is set by the ratio of the Doppler shift to the subcarrier spacing: fc/Δf is smaller at 28 GHz / 120 kHz than at 5 GHz / 15 kHz, so 300 km/h is still below the 16-QAM threshold of 13.1%, while 5 GHz / 15 kHz exceeds it at about 230 km/h.

#### E2 [zh] E2 樣本太少
- look: 快照數 L=4、天線數 N=8，R̂ 的秩只有 4，診斷區會出現「秩虧」。SMI（偽逆）平均 SINR 約 22.46 dB，DL 約 27.81 dB，BEAMSPACE 約 26.36 dB。
- why: 樣本太少時協方差矩陣不能求逆：SMI 改用偽逆（只用有資料的方向），DL 在對角線加一點量讓它可逆，BEAMSPACE 先把 8 根天線壓成 3 個波束再求逆。
- real: 通道變化很快時能收集的快照有限，所以小樣本是常態；在這個設定下 DL 比 SMI 高 5.3 dB，BEAMSPACE 高 3.9 dB。

#### E2 [en] E2 Too few snapshots
- look: With L=4 snapshots and N=8 antennas the rank of R̂ is only 4 and the diagnosis shows "rank deficient". The mean SINR is about 22.46 dB for SMI (pseudo-inverse), 27.81 dB for DL and 26.36 dB for BEAMSPACE.
- why: With too few samples the covariance matrix cannot be inverted: SMI uses the pseudo-inverse (only the directions that have data), DL adds a little to the diagonal so that it becomes invertible, BEAMSPACE first squeezes the 8 antennas into 3 beams and then inverts.
- real: When the channel changes quickly only a few snapshots can be collected, so small samples are the normal case; in this setting DL is 5.3 dB and BEAMSPACE 3.9 dB above SMI.

#### E3 [zh] E3 自我抵消
- look: 訓練資料含訊號（MPDR）、指向偏差 3°：SMI 的 SINR 從 -5.56 dB（0°）掉到 -19.75 dB，DL 從 13.98 dB 掉到 -0.27 dB，MMSE-P 維持 28.58 dB。對照：訓練資料不含訊號時，SMI 在 0° 有 28.57 dB。
- why: MPDR 要在「不扭曲指定方向」的限制下讓輸出功率最小；訓練資料裡有目標時，只要指定的方向有一點偏差，演算法就把真正的目標當成干擾壓掉（自我抵消）。MMSE-P 用已知導頻直接量出目標，不用你指定的方向，所以不受影響。
- real: 接收機通常只有含訊號的訓練資料，而角度誤差（天線校正、DOA 估計、目標移動）一定存在。DL 整體比 SMI 高約 19.5 dB，但兩者都下降了約 14.2 dB，真正不受影響的是用導頻的 MMSE-P。

#### E3 [en] E3 Self-nulling
- look: Training data with the signal (MPDR) and a 3° pointing error: the SINR of SMI falls from -5.56 dB (0°) to -19.75 dB, DL from 13.98 dB to -0.27 dB, while MMSE-P stays at 28.58 dB. For comparison, without the signal in the training data SMI has 28.57 dB at 0°.
- why: MPDR minimises the output power under the constraint of not distorting the assumed direction; when the target is in the training data and the assumed direction is slightly off, the algorithm suppresses the real target as if it were interference (self-nulling). MMSE-P measures the target from the known pilots and does not use the assumed direction, so it is not affected.
- real: A receiver usually has only training data that contain the signal, and angle errors (antenna calibration, DOA estimation, target motion) are unavoidable. DL is about 19.5 dB above SMI, but both fall by about 14.2 dB; only MMSE-P, which uses pilots, is really unaffected.

#### E4 [zh] E4 指向偏差
- look: 圖中是 SMI 的 SINR 對指向偏差 δθ。訓練資料不含訊號時緩慢下降：0° 28.59 dB、3° 27.96 dB、5° 26.8 dB；含訊號時急遽下降：0° -5.72 dB、1° -12.03 dB、3° -19.89 dB。
- why: 不含訊號時，偏差只是讓波束稍微偏離目標，增益慢慢變小；含訊號時，偏差一出現就觸發自我抵消，所以只差 1° 也會掉很多。
- real: 你指定的方向通常來自 DOA 估計（面板 F）或天線校正，兩者都有誤差，這個誤差就是這裡的 δθ。

#### E4 [en] E4 Pointing error
- look: The plot shows the SINR of SMI against the pointing error δθ. Without the signal in the training data it falls slowly: 28.59 dB at 0°, 27.96 dB at 3°, 26.8 dB at 5°; with the signal it falls steeply: -5.72 dB at 0°, -12.03 dB at 1°, -19.89 dB at 3°.
- why: Without the signal the error merely points the beam slightly away from the target, so the gain shrinks slowly; with the signal the error triggers self-nulling at once, so even 1° costs a lot.
- real: The direction you assume usually comes from a DOA estimate (panel F) or from antenna calibration; both have errors, and that error is the δθ here.

#### E5 [zh] E5 通道老化
- look: 更新延遲 τ 從 0 增加到 10 ms（車速 300 km/h、離軌道 30 m）：SMI 的 SINR 從 28.31 dB 降到 25.48 dB，DL 從 28.31 dB 降到 23.75 dB。
- why: 權重是 τ 之前估計的；這段時間車子往前開，干擾源相對於車的角度改變了，零陷還停在舊角度，沒有對準干擾。
- real: 從量測到套用權重一定有處理延遲；離軌道越近，角度變得越快（離軌道 5 m 時 DL 的 SINR 從 26.19 dB 降到 17.17 dB），所以高速又靠近時要更常更新權重。

#### E5 [en] E5 Channel aging
- look: Increase the update latency τ from 0 to 10 ms (300 km/h, 30 m from the track): the SINR of SMI falls from 28.31 dB to 25.48 dB, that of DL from 28.31 dB to 23.75 dB.
- why: The weights were estimated τ earlier; meanwhile the vehicle has moved on, the angle of the jammer as seen from the vehicle has changed, and the null still points at the old angle instead of at the jammer.
- real: There is always a processing delay between measuring and applying the weights; the closer the track, the faster the angle changes (at 5 m DL falls from 26.19 dB to 17.17 dB), so at high speed and short distance the weights must be updated more often.

#### E6 [zh] E6 速度如何變成指向偏差
- look: 車速從 0 掃到 300 km/h（離軌道 10 m、τ = 10 ms、目標 20°、干擾 −30°）：有效指向偏差 δθ_eff 增加到 0.54°（只有 3 dB 主瓣寬 14° 的 3.89%），干擾的漂移是 -1.15°。SMI（不含訊號）的 SINR 從 28.61 dB 降到 27.36 dB，DL 降到 23.91 dB，MMSE-P 從 28.63 dB 掉到 -11.17 dB。含訊號訓練的 SMI 趨勢不同：從 -5.61 dB 先升到 7.14 dB（175 km/h），再降到 5.05 dB。同樣的速度下，頻域的 ICI 底限 N_ICI/S 是 -16.08 dB（100 km/h：-25.58 dB），比 SMI（不含訊號）因老化多出來的損失（300 km/h：-33.38 dB）大；掃描的第一個非零速度（25 km/h）起 ICI 就已經比較大，在 0–300 km/h 內沒有交叉點。以上是 fc = 5 GHz、Δf = 15 kHz、離軌道 10 m、τ = 10 ms 的結果；在 fc/Δf ∈ {5 GHz/15 kHz, 28 GHz/120 kHz}、離軌道 ∈ {10, 30, 100, 500} m、τ ∈ {1, 2, 10} ms、v ∈ {25, 50, 100, 200, 300} km/h、SMI 與 DL 共 240 格中，老化損失 ≥ ICI 損失的有 0 格；最接近的是 28 GHz/120 kHz、10 m、10 ms、300 km/h、SMI（老化 5.53 dB，ICI 9.91 dB）。
- why: EVM 上升來自 ICI（頻域：都卜勒讓子載波失去正交，看 E1）；SINR 下降來自老化（空間域：權重是 τ 之前估計的，目標與干擾的角度已經移動，通道相位也轉過了）。把 τ 設成 0（沒有角度偏差）時，SMI（不含訊號）在 300 km/h 仍有 28.51 dB，所以它的損失來自 τ 這段延遲。含訊號訓練的 SINR 在 0 km/h 很低，是因為權重把訓練資料裡的目標當成干擾壓掉（自我抵消）；速度提高後通道在訓練窗內轉動，自我抵消變弱（τ = 0、300 km/h 時 6.58 dB，100 km/h 時 0.11 dB），所以 SINR 上升，並不是指向變準了。MMSE-P 靠導頻估計互相關，通道一旋轉就失效。
- real: 要分清楚是哪一種損失：EVM 的底限看頻域（換較大的子載波間隔或較低載頻），SINR 的損失看空間域（縮短更新延遲、離軌道遠一點，或用對角加載的方法）。指向模式選「由移動性決定」時，δθ 滑桿失效，偏差由 d_min、v、τ 自動決定。

#### E6 [en] E6 How speed becomes a pointing deviation
- look: Sweep the speed from 0 to 300 km/h (10 m from the track, τ = 10 ms, target 20°, jammer −30°): the effective pointing deviation δθ_eff grows to 0.54° (only 3.89% of the 14° 3 dB beamwidth), the jammer drifts by -1.15°. The SINR of SMI (training without the signal) falls from 28.61 dB to 27.36 dB, DL to 23.91 dB, MMSE-P drops from 28.63 dB to -11.17 dB. SMI trained with the signal behaves differently: it first rises from -5.61 dB to 7.14 dB (at 175 km/h) and then falls to 5.05 dB. At the same speeds the frequency-domain ICI floor N_ICI/S is -16.08 dB (100 km/h: -25.58 dB), larger than the extra loss of SMI (no signal) caused by aging (300 km/h: -33.38 dB); the ICI floor is already the larger one at the first non-zero speed of the sweep (25 km/h), so there is no crossing between 0 and 300 km/h. These are the results for fc = 5 GHz, Δf = 15 kHz, 10 m from the track and τ = 10 ms; over the grid fc/Δf ∈ {5 GHz/15 kHz, 28 GHz/120 kHz}, distance ∈ {10, 30, 100, 500} m, τ ∈ {1, 2, 10} ms, v ∈ {25, 50, 100, 200, 300} km/h, SMI and DL (240 cells) the aging loss reaches the ICI loss in 0 cells; the closest is 28 GHz/120 kHz, 10 m, 10 ms, 300 km/h, SMI (aging 5.53 dB, ICI 9.91 dB).
- why: The EVM rises because of ICI (frequency domain: the Doppler shift destroys the orthogonality of the subcarriers, see E1); the SINR falls because of aging (spatial domain: the weights were estimated τ earlier, the target and the jammer have moved and the channel phase has rotated). With τ = 0 (no angle deviation) SMI (no signal) still reaches 28.51 dB at 300 km/h, so its loss comes from the delay τ. The SINR with the signal in the training data is low at 0 km/h because the weights suppress the target in the training data as if it were interference (self-nulling); at higher speed the channel rotates inside the training window and the self-nulling weakens (τ = 0: 6.58 dB at 300 km/h, 0.11 dB at 100 km/h), so the SINR rises, not because the pointing gets better. MMSE-P estimates the cross-correlation from pilots and fails as soon as the channel rotates.
- real: Tell the two losses apart: the EVM floor is a frequency-domain matter (larger subcarrier spacing or lower carrier frequency), the SINR loss is a spatial-domain matter (shorter update latency, a larger distance to the track, or diagonal loading). With the pointing mode set to "set by mobility" the δθ slider is inactive and the deviation follows from d_min, v and τ.

## 7. 瀏覽器驗證

- 方法：Edge＋playwright-core，繁中與英文各一輪（預設、舊模型、withSignal、28 GHz/120 kHz、指向 3°、DOA 兩種訓練、掃描、E0–E6、1280 與 1024 寬度）；檢查截斷／重疊文字、英文模式殘留中文、空白畫布、console 錯誤、CSV 欄位。最終結果：**0 個問題**。
- 已修問題：(1) 英文模式下殘留的未翻譯中文片段被重複改寫，使 MutationObserver 無窮迴圈、頁面卡死（寫入前比對、監聽 title、補短語表，T19a 擴充）；(2) 標籤截斷與重疊（`.kv .k`、`.p-head`、`.lab` 換行）；(3) 1024 寬度下 D 區標籤被硬拆行如「(geni e)」（`overflow-wrap: break-word`、欄寬 190 px）。
- 尚存的小觀察：頁首管線列的 v、f_d 讀數最多落後 250 ms 更新，截圖若剛好在載入後瞬間拍攝會看到舊值；不影響物理與 CSV。
- 截圖（`docs/screenshots/`，每項 `-zh`、`-en` 各一，共 37 個檔案含 `report.json`）：`00-default`、`01-legacy`、`02-withsignal`、`03-28ghz-120khz`、`04-pointing-3deg`、`05-doa-withsignal`、`06-doa-signalfree`、`07-sweep`、`08-width-1280`、`08-width-1024`、`exp-e0`、`exp-e1-28`、`exp-e1-5`、`exp-e2`…`exp-e6`。與上一版相比，37 個檔案內容都變了（版面調整改變欄寬，且畫面含動畫時間）。

## 8. MATLAB 部分（**未在本環境執行**）

未執行的腳本：`matlab/selftest.m`、`check_my_work.m`、`verify_cases.m`、`mc_independent.m`、`plot_figures.m`、`lab_mmse_mvdr_demo.m`、`exercises/ex1–ex5`、`solutions/ex*_solution.m`。Node 端只驗證了：案例資料可由匯出資料重算（T20b，權重逐位相同）、批改邏輯對參考解答 5/5、對常見錯誤判為未通過（T20f）。

使用者執行順序（目前資料夾切到 `matlab/`，只需基本 MATLAB、R2016b 以上）：
1. `selftest` → 預期「通過 5 / 5 關」。
2. 依序做 `exercises/ex1…ex5`（先自己寫，卡住再看提示，最後看解答），每做完一關執行 `check_my_work`，最後一行「通過 X / 5 關」。
3. `verify_cases` → 23 案例 × 6 演算法＋E0 的 9 項檢查，共 147 項，每行 PASS/FAIL；L=N 的案例可能病態（κ 到 1e8），誤差與估計值同量級時是數值條件問題。
4. `mc_independent` → 40 列表，最後一行 `m = 40，Bonferroni 臨界值 z_crit …  ->  PASS`。
5. `plot_figures` → `figs/` 出現 fig1–fig4（pdf＋png）。
6. `lab_mmse_mvdr_demo`（可貼在講義程式後面）→ 理論協方差的 MMSE 與 MVDR A、B 的 SINR 與原始 EVM。
若任何一步報錯，請把畫面貼回來。

## 9. 已知限制與不確定處

- **MMSE-P 無頻偏補償**：v ≥ 30 km/h（視窗內 LoS 旋轉約 1 圈以上）SINR 崩潰（θ₁=0°、K=20 dB：v=0 為 28.6 dB，30 km/h −20.0 dB，100 km/h −9.7 dB，300 km/h −12.4 dB）。原因是 r̂_xd = (1/L)Σ x conj(s) 對公共相位旋轉做相干平均（genie 去旋轉可回到 28.5 dB）。規格的 lag-1 估計器在 SIR −10 dB 下不可用；週期圖法未通過 T15a/T15b（第 3 項）。
- **T13a 的 ICI 高斯近似**：符元層的 ICI 項以高斯雜訊生成（功率取自同一個 N_ICI/S），T13a 只驗證鏈路的算術（縮放、正規化、樣本估計量），不驗證 ICI 物理（那是 T4／T5）。signalFree、v=300 那一格在單一種子區塊曾達 z=3.72，其餘區塊與 2000 次平均都沒有偏差，原因未能進一步指認。
- **積分相位建模**：直線軌道、精確幾何 θ̇ = v sinθ|sinθ|/d_min，積分都卜勒相位 φ = 2π[R(0)−R(t)]/λ；θ₁ = 0（在軌道上）不漂移；不含轉彎、加速、多普勒的二階以外效應。
- **干擾源單一路徑**：干擾只有一條直達路徑（無擴散）；干擾波形預設高斯（OFDM 近似），`qpsk` 為舊行為；jammer 波形選擇改變亂數序列與部分數值（golden Commit 19）。
- **權重估計時刻**：權重用 t_app−τ 的角度，漂移從未被補償；`pointingMode = 'mobility'` 與 δθ=0 逐位相同（不是新物理）；`manual` 時總偏差為 δθ_eff − δθ。
- **E6 的峰值位置很平**：withSignal 的 SMI 在 100–150 km/h 都落在 6.7–7.2 dB；文字寫的 175 km/h（7.14 dB）是固定種子下的最大值，S3 用不同種子最高點在 125 km/h。
- **E6 的適用範圍**：S1 網格（fc/Δf ∈ {5 GHz/15 kHz, 28 GHz/120 kHz}、d_min ∈ {10, 30, 100, 500} m、τ ∈ {1, 2, 10} ms、v ∈ {25…300} km/h、SMI／DL）240 格中老化損失都小於 ICI 損失，最接近的是 28 GHz/120 kHz、10 m、10 ms、300 km/h、SMI（老化 5.53 dB、ICI 9.91 dB）；MMSE-P 與 withSignal 不在該網格內。
- **L<N 時 MMSE-P 的 MATLAB 對照（`verify_cases.m`，案例 `static_L4_signalFree_d3`）**：L<N 時 MMSE-P 以偽逆求解，網頁與 MATLAB 的權重差異來源為捨入敏感度（400 個隨機 L=4 實現中 14% 單靠 R̂ 加總順序改變即超過 1e-9）；該案例的 FAIL 為已知、已調查、不影響任何結論。兩邊的偽逆閾值相同（1e-10·λmax），保留部分的條件數約 1e4。`verify_cases.m` 現在把 L<N 的案例標為 `rank-deficient (L<N)` 並在總結中分開計數，判準數值（1e-9）未改，FAIL 仍顯示為 FAIL。調查腳本與結果在 `docs/diagnostics/l4_*`，MATLAB 端未在本環境執行。
- **其他**：L=N 附近 R̂ 病態（MMSE-P、SMI 都有低谷，是估計量性質而非求逆實作，B1）；MUSIC 取最接近名義角的峰，網格 0.5° 量化；E0 使用 K = 400 dB（超出滑桿範圍 20 dB）；原始 EVM 的 LoS 公共相位已扣除（講義通道沒有隨機載波相位）；講義公式編號未知，只對照變數名稱與慣例（θ_lab = θ_sim）；預設模型與訓練模式只在 Commit 17 切換。
- **關於 S1–S4 敏感度調查**：S1–S4 調查已於本報告撰寫前完成並 commit（結果在 `docs/diagnostics/`），本報告僅引用既有結果，未再執行新的調查。

## 10. 可用於論文／報告的結論句（僅列有數據支持的；括號為適用範圍與依據）

設定縮寫：N = 8、d = 0.5λ、SNR 20 dB、SIR −10 dB、θ₂ = 40°（除非另註）、unified 單一實現模型、Gaussian 干擾。

核對說明：本項每個數字已由 `tests/check_report_numbers.js` 對照其來源核對，結果在 `docs/diagnostics/report_number_check.txt`。其中 T17b2、T17c、T17d 的數字（第 1、2、10 句，以及第 2 句的範圍與差值）來自**已存檔的完整套件輸出** `docs/diagnostics/suite_output_2026-10-06_198of198.txt`（2026-10-06，198/198 那次），這次沒有重新執行 T17；T11b 與 T18c 的數字來自新鮮重跑（`docs/diagnostics/fresh_T11_output.txt`、`fresh_T18_output.txt`）。

1. 在 fc = 5 GHz、Δf = 15 kHz、d_min = 30 m、L = 100、K = 20 dB、v = 0、τ = 0、θ₁ = 0°（θ₂ = 40°、SNR 20 dB、SIR −10 dB、N = 8）時，不含訊號訓練的 SMI 與 MMSE-P 的 SINR 都約 28.6 dB，距遺傳界 SINR_opt 約 0.3 dB。（MMSE-P：T11b，2000 次實現；SMI：E4 的 0° 點 28.59 dB 與 T17c，1000 次實現；T14 為固定種子的回歸比對，不是這些平均值的來源；v = 0 時 d_min 不影響結果。）
2. 同一設定下（fc = 5 GHz、Δf = 15 kHz、v = 0、τ = 0、θ₁ = 0°、θ₂ = 40°、L = 100、K = 20 dB；d_min 在 v = 0 時不影響）若訓練資料含目標（MPDR，withSignal）、δθ = 0，SMI 的 SINR 降到約 −5.5 至 −5.7 dB（自我抵消）；在 withSignal 下指向偏差 3° 時 SMI 與 DL 比 0° 低 10 dB 以上，MMSE-P 不受指向偏差影響。（T17d −5.5、E3 −5.56、E4 −5.72（T17d 的輸出只有一位小數），三者均為 1000 次實現、種子不同；各平均值的標準誤差約 0.10–0.11 dB，兩兩最大差 0.22 dB，約為差的標準誤差（0.15 dB）的 1.5 倍，彼此在統計上相容（`docs/diagnostics/withsignal_se_check.txt`）；signalFree 時 3° 只降約 0.6 dB（SMI 28.6→28.0，T17d），不屬於本句。）
3. 以講義的理論協方差（N=8、θ₁=−20°、θ₂=30°、SNR 30 dB、SIR 0 dB、K=∞（400 dB）、v = 0、τ = 0、QPSK；fc = 5 GHz、Δf = 15 kHz、d_min = 30 m 為預設值，在 v = 0 時不影響），MMSE-M、MVDR 方法 A、B 的 SINR 都等於 SINR_opt（38.97 dB）；改用 L=1000 的樣本協方差且訓練含訊號時降為 21.88 dB，不含訊號的 SMI 仍有 38.94 dB。（T21，預設 300 次實現；E0 數字為 200 次實現的平均；L = 1000。）
4. 快照數 L=4 < N=8 時秩只有 4，DL（27.81 dB）與 BEAMSPACE（26.36 dB）高於 SMI（22.46 dB）。（E2，1000 次實現；fc = 5 GHz、Δf = 15 kHz、v = 0、τ = 0、θ₁ = 0°、θ₂ = 40°、K = 20 dB、SNR 20 dB、SIR −10 dB、signalFree；DL 為 γ_rel = +10 dB；d_min 在 v = 0 時不影響。）
5. 在 v = 300 km/h、θ₁ = 0°、θ₂ = 40°、d_min = 30 m、τ = 0、L = 100、K = 20 dB（SNR 20 dB、SIR −10 dB、signalFree）下，28 GHz/120 kHz 的 SMI EVM 為 12.31 %、5 GHz/15 kHz 為 17.12 %；後者在約 230 km/h 越過 16-QAM 的 13.14 % 門檻，前者在 0–300 km/h 內沒有越過。（E1，SMI 的 EVM 取各實現均方根，每點 300 次實現、速度步長 10 km/h；越過速度的解析度為 10 km/h。）
6. 更新延遲 τ 從 0 增到 10 ms（v=300 km/h、d_min=30 m、θ₁ = 0°、θ₂=−30°、K=20 dB、L=100、fc = 5 GHz、Δf = 15 kHz、SNR 20 dB、SIR −10 dB、signalFree）使 SMI 的 SINR 從 28.31 降到 25.48 dB、DL 降到 23.75 dB；離軌道 5 m 時 DL 從 26.19 降到 17.17 dB。（E5，各 1000 次實現；T3a 在 K→∞、θ₁ = 0°、θ₂ = −30°、d_min = 5 與 30 m 下只驗證「SINR 隨 τ 單調下降」的方向，不是這些數字的來源。）
7. 在 5 GHz/15 kHz、d_min ≤ 500 m、τ ≤ 10 ms、v ≤ 300 km/h 與 28 GHz/120 kHz（Δf = 120 kHz）的對應網格內（d_min ∈ {10, 30, 100, 500} m、τ ∈ {1, 2, 10} ms、v ∈ {25, 50, 100, 200, 300} km/h；SMI、DL（γ_rel = +10 dB），signalFree，K=20 dB、θ₁=20°、θ₂=−30°、L=100、N=8、SNR 20 dB、SIR −10 dB），空間域的老化損失都小於 ICI 損失；最接近的情況是 28 GHz/120 kHz、10 m、10 ms、300 km/h：老化 5.53 dB 對 ICI 9.91 dB。（S1，240 格、每格 500 次實現：`docs/diagnostics/s1_sensitivity.csv`、E6；不含 MMSE-P 與 withSignal，網格以外不外推。）
8. 含訊號訓練的 SMI 的 SINR 在高速下先升後降，原因是擴散路徑間的都卜勒差：K→∞ 時不上升（v=0 為 11.4 dB，300 km/h 為 6.7 dB），K=20 dB 時從 −5.5 dB 升到 7.2 dB（125 km/h）後降到 4.6 dB（300 km/h）。（適用範圍為 E6 設定：fc = 5 GHz、Δf = 15 kHz、d_min = 10 m、τ = 10 ms、θ₁ = 20°、θ₂ = −30°、L = 100、N = 8、SNR 20 dB、SIR −10 dB、withSignal 的 SMI；K→∞（400 dB）與 K = 20 dB 的對照取自 S3，每點 500 次實現（目標協方差特徵值佔比為 200 次），固定參數即上述 E6 設定、速度 0–300 km/h 步長 25 km/h；結論是「K→∞ 時不上升」，不含 MMSE-P 或其他 d_min、τ。對應 E6、T22c（τ = 0 對照）。）
9. 在 θ₁ = 0°、K = 20 dB、v = 30–300 km/h、L = 100（fc = 5 GHz、Δf = 15 kHz、d_min = 30 m、τ = 0、θ₂ = 40°、SNR 20 dB、SIR −10 dB）下，沒有頻偏補償的 MMSE-P 在 30 km/h 起崩潰（−20.0 dB），用真實 LoS 相位去旋轉（genie）可回到 28.5 dB（300 km/h；30 km/h 為 28.6 dB）。θ₁ = 45° 時，用真實相位去旋轉在 300 km/h 只回到 14.76 dB（K = 20 dB），此結論不能推廣到其他角度。（B4，每格 500 次實現：`docs/diagnostics/b4_afc.txt`；MMSE-P 的 R̂ 含目標，signalFree 與 withSignal 相同；這是診斷不是演算法結果。）
10. 樣本 SMI（L = 20000，signalFree，δθ = 0° 與 1°）與母體解的差，與 Reed–Mallett–Brennan 期望損失 −0.0015 dB 相符（z = −2.07、−1.78）。（T17b、T17b2，各格 100 次實現；N = 8、K = 20 dB、v = 0、τ = 0、θ₁ = 0°、θ₂ = 40°、SNR 20 dB、SIR −10 dB、QPSK；fc = 5 GHz、Δf = 15 kHz、d_min = 30 m 為預設值，v = 0 時不影響；δθ ≥ 3° 公式不嚴格適用。）
11. MUSIC 與 Capon 的目標角估計誤差在 L ≥ 12（θ₁ = 0°、θ₂ = 40°）與源分離 ≥ 20°（θ₁ = 0°、L = 100）時約 0.2°（受 0.5° 網格限制）。（T18c，各格 500 次實現；N = 8、K = 20 dB、SNR 20 dB、SIR −10 dB、v = 0、τ = 0、withSignal、雙源模型；fc = 5 GHz、Δf = 15 kHz、d_min = 30 m 為預設值，v = 0 時不影響。）

不列入論文結論的：T13a（鏈路算術檢查）、MATLAB 部分（未執行）、withSignal 與 MMSE-P 的高速行為（僅為診斷）。
