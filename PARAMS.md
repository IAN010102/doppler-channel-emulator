# PARAMS.md — parameters, conventions and the two channel models

`core.js` holds the DOM-free physics (browser global `Core`, or `require('./core.js')` in Node); `index.html` is UI + drawing only.
All defaults live in `CONFIG` (top of `core.js`); all randomness goes through `Core.rng()` (mulberry32, `Core.setSeed(s)`).
The UI picks a fresh random seed at every page load.

## 1. CONFIG defaults

| key | value | meaning |
|---|---|---|
| `fc` | 5 GHz | carrier frequency (UI: 3.5 / 5 / 7 / 15 / 28 GHz or custom; §8) |
| `c` | 3×10⁸ m/s | speed of light |
| `scs` | 15 kHz | OFDM subcarrier spacing Δf (UI: 15 / 30 / 60 / 120 kHz; §8) |
| `N` | 8 | default number of array elements (UI: 2–16) |
| `d_lambda` | 0.5 | element spacing in wavelengths (ULA) |
| `cpRatio` | 0.07 | CP ratio; snapshot period `T_snap = (1 + cpRatio)/Δf` = 71.33 µs (5G-NR-like symbol length) |
| `d_min` | 30 m | perpendicular distance from the ground point to the straight track (§4); renamed from `R_min` |
| `SIGMA_ANG_DEG` | 10° | angular spread σ_θ of the diffuse paths |
| `M_SCAT` | 8 | diffuse paths per snapshot — **legacy model only** |
| `M_UNIFIED` | 32 | diffuse paths per trial — **unified model** |
| `LAMBDA_Q`, `REL_Q` | 30, 10 | loading rule of the tapered (GSC) adaptive path (unchanged) |
| `REFRESH` | 0.2 | sliding-window refresh fraction — **legacy model only** |
| `iciWarnDb`, `iciSevereDb` | −30, −20 dB | N_ICI/S thresholds of the ICI diagnosis (§8) |
| `gammaRelDb` | +10 dB | unified model, DL: γ = 10^(γ_rel/10)·σ_n² (see §7); legacy keeps absolute γ = 0.01 |
| `trainMode` | `'withSignal'` | `'withSignal'` (MPDR) or `'signalFree'` (MVDR); see §6; page flag `?train=signalFree` |
| `model` | `'legacy'` | `'legacy'` \| `'unified'`; the page also accepts `?model=unified` |

Other defaults (UI): L = 100, θ₁ = 0°, θ₂ = 40°, SNR = 20 dB (per element), SIR = −10 dB, v = 0, K = 20 dB, σ_φ = 0°, τ = 0, γ = 0.01 (absolute), 16-QAM.
Not modelled (unchanged): number of OFDM subcarriers, sampling rate, time-domain waveform, FFT size.

## 2. Angle and Doppler conventions (unchanged from the existing code)

- θ is measured from the array broadside (the x axis). The array lies along y; the receiver/target heading is **+x**, so θ is also the angle between a path and the direction of travel.
- Steering vector `a_n(θ) = exp(−j2π n d sinθ)`, n = 0…N−1. The true array is `Γ ⊙ a(θ)` with `Γ = diag(e^{jφ_n})` (phase mismatch); weights always use the nominal `a`.
- `f_m = v·fc/c`, path Doppler `f_d,i = f_m·cosθ_i`, `ε_i = f_d,i/Δf`, `sinc(x) = sin(πx)/(πx)`.
- Diffuse angle distribution (existing code): **Gaussian**, θ_m = θ₁ + δ_m, δ_m ~ N(0, σ_θ²), σ_θ = 10°. Complex gains: CN(0, 1/((K+1)M)) per path (legacy: re-drawn every snapshot; unified: drawn once per trial).

## 3. Unified model (`model = 'unified'`)

Path set per trial (one draw of {φ₀, δ_m, g_m}, parameter-free; K, v, τ, θ₁ are applied when the draw is used):

| path | angle | complex gain β |
|---|---|---|
| LoS | θ₁ | √(K/(K+1))·e^{jφ₀}, φ₀ ~ U[0, 2π) |
| diffuse m = 1…M | θ₁ + δ_m | g_m, g_m ~ CN(0, 1/((K+1)M)) |
| jammer | θ₂ | √P_j (single LoS path; unit-power QPSK j_n) |

Time: snapshot n at `t_n = n·T_snap` (n = 0…L−1); estimation time `t_est = (L−1)·T_snap`; application time `t_app = t_est + τ` (τ = "update latency", same meaning as before: weights computed from the window, applied τ later).
Reference angles are the slider values = angles **at t_app**; the path angle follows the exact straight-track law θ_i(t) = trackAngle(θ_i,0, v, d_min, t − t_app), θ̇ = v sinθ|sinθ|/d_min for target paths and jammer alike (§4). The Doppler phase is the integral along the track, φ_i(t) = 2π[R_i(0) − R_i(t)]/λ (§4); the ICI uses the instantaneous `f_d,i = f_m cosθ_i,0` at t_app.

Channel and snapshots:

```
h(t)  = Σ_i β_i · Γ⊙a(θ_i(t)) · exp(jφ_i(t))                (φ_i: integrated Doppler phase, §4)
x_n   = h(t_n)·s_n + √P_j · Γ⊙a(θ₂(t_n))·exp(jφ₂(t_n))·j_n + noise_n        (s_n, j_n unit-power QPSK, noise CN(0,σ²))
```

The covariance estimate and all five weight designs (FOURIER / MMSE / SMI / DL / BEAMSPACE, taper variants) are unchanged and use `R̂ = (1/L)Σ x_n x_nᴴ`, with the look direction `a(θ̂₁)`, θ̂₁ = θ₁(t_est).

Per-trial metrics at `t_app`:

```
SINR_inst = |wᴴ h(t_app)|² / ( P_j |wᴴ a(θ₂)|² + σ² ‖w‖² )
q_i       = |wᴴ a(θ_i)·β_i|² / Σ_k |wᴴ a(θ_k)·β_k|²
N_ICI/S   = Σ_i q_i (1 − sinc²(ε_i))
EVM       = √( 1/SINR_inst + N_ICI/S )
```

**Approximations / assumptions (stated, not hidden):**
1. The per-path leakage terms are added *non-coherently* (cross terms between paths are dropped); `1 − sinc²(ε)` is the N_FFT → ∞ limit of `1 − |sin(πε)/(N sin(πε/N))|²`.
2. The receiver performs **no Doppler compensation** (no AFC): every path keeps its own `ε_i`.
3. ICI is treated as additional Gaussian noise relative to the output signal power; SER still uses the existing Gaussian/M-QAM closed form on that EVM.
4. The ICI floor is `EVM_floor = √(Σ q_i (1 − sinc²(ε_i)))` with `q_i = |β_i|²/Σ|β_k|²` (no spatial filtering), evaluated as the **expectation over the path distribution** in closed form/numerically: `K/(K+1)(1 − sinc²(ε_LoS)) + 1/(K+1)·E_δ[1 − sinc²(ε_m cos(θ₁+δ))]` (ratio of expected powers).

UI semantics of a "trial": the path realisation persists until RESET (or a change of `M_UNIFIED`); symbols and noise are fresh at every update. Tests and the benchmark set `freshRealization = true`, i.e. one realisation per call.

Deliberate differences from the legacy model (spec-driven): static fading within a trial (no per-snapshot re-draw); the jammer is a single path (legacy: Rician with the shared K); the metric is instantaneous `SINR_inst`, not an expected-covariance SINR.

## 4. Geometry: angle drift, `d_min` and the Doppler sign (unified model)

**Setting.** A fixed ground point (target or jammer) lies at perpendicular distance `d_min` from a straight track. The receiver moves along **+x** at speed v. θ is the angle between the line of sight (receiver → source) and the heading +x (the existing convention: θ measured from the array broadside = heading, `f_d = f_m cosθ`, `a_n = e^{−j2π n d sinθ}`; θ > 0 on one side of the track, θ < 0 on the other; |θ| < 90° = source ahead).

**Derivation.** Put the source at (x_s, d_s), the receiver at (x_r(t), 0), x_r = vt, r = √((x_s−x_r)² + d_s²), d_min = |d_s|:

- cosθ = (x_s − x_r)/r, sinθ = d_s/r. Define c = cosθ/|sinθ| = (x_s − x_r)/d_min (valid on both sides of the track).
- dc/dt = −v/d_min (the receiver moves v·dt along x, the perpendicular distance does not change). For θ > 0, c = cotθ.
- Since dc/dt = −θ̇/sin²θ·sign(sinθ):  **θ̇ = v·sinθ·|sinθ| / d_min** (= +v sin²θ/d_min for θ > 0). Range rate ṙ = −v cosθ (approaching for θ < 90°).
- Doppler: `f_d = f_m cosθ`, `ḟ_d = −f_m sinθ·θ̇ = −f_m v sin²θ|sinθ|/d_min ≤ 0` for θ > 0: as the receiver closes on the point, θ grows and f_d falls (positive f_d = approaching, consistent with the sign convention f_d = +f_m cosθ). So the drift of θ and the sign of f_d are consistent in this convention.
- Exact solution used in the code (`Core.trackAngle`): `θ(t + Δt) = sign(sinθ)·atan2(1, cosθ/|sinθ| − vΔt/d_min)`. A source exactly on the track (θ = 0 or 180°) does not drift. Checked against explicit receiver/source coordinates by `tests/t7_geometry.js`.

**Difference from the previous code.** The previous unified drift was `θ̇ = −v sinθ/R_min` (jammer `+v sinθ/R_min`): wrong magnitude (equal to the exact value only at |θ| = 90°, over-estimated by 1/|sinθ| elsewhere) **and opposite sign** for the target in this angle convention (θ moves away from the heading as the receiver approaches). The variable is now named `d_min` (`CONFIG.d_min`, `Sys.d_min`, default 30 m).

**Target and jammer** use the same law, each with its own current angle (θ_i at t_app for the path table; `thTo`/`thJo`, the look and null angles at t_est, are `trackAngle(θ, v, d_min, −τ)`). The "jammer moves along −x" special case is removed (a fixed ground jammer seen from a receiver moving along +x has the same geometry as the target).

**Legacy model**: kept as it was (`thTo = θ₁ + v sinθ₁ τ/d_min`, `thJo = θ₂ − v sinθ₂ τ/d_min`, only the variable was renamed from `R_min`).

**Integrated Doppler phase (Commit 9; replaces the first-order phase 2π·f_d(t_app)·t).** For every path (LoS, each diffuse path, jammer):

```
φ_i(t) = 2π [R_i(0) − R_i(t)] / λ,     R_i(t) = d_min / |sin θ_i(t)|,     θ_i(t) = trackAngle(θ_i,0, v, d_min, t − t_app)
```

`θ_i,0` is the angle at t_app (the slider angle for the target LoS and the jammer; θ₁ + δ_m for diffuse paths). t = 0 is the first training snapshot, so φ_i(0) = 0 and the estimation window starts at t = 0 as before. Closed form (stable for d_min → ∞): with c = cosθ/|sinθ| − v(t − t_app)/d_min, `R(ta) − R(tb) = v (tb − ta)(c_a + c_b)/(√(1 + c_a²) + √(1 + c_b²))` (`Core.trackPhase`). `dφ/dt = 2π f_m cosθ(t) > 0` while approaching, and φ → 2π f_d t for d_min → ∞ (relative deviation ≈ v·t/d_min: 1.6×10⁻⁶ at 10⁶ m). Verified in `tests/t8_doppler_phase.js`.

**Modelling assumption (stated, not derived):** every path — LoS, diffuse and jammer — is treated as a **far-field plane wave that shares the same d_min**; a diffuse path is the same ground-point geometry with its own initial angle θ_i,0 = θ₁ + δ_m. (In reality a scatterer has its own distance; this is not modelled.) A source exactly on the track (sinθ = 0) uses R(ta) − R(tb) = v(tb − ta)cosθ.

**ICI.** `ε_i = f_d,i/Δf` still uses the instantaneous `f_d,i = f_m cosθ_i,0` **at t_app** (the application time): the angle change within one OFDM symbol (T_snap ≈ 71 µs) is negligible (≈ 0.006° at 300 km/h, d_min = 30 m), so the symbol-level ICI is evaluated with the Doppler at that instant, while the slow phase evolution over the training window and τ is carried by φ_i(t).

**Diagnostic (panel D, `Sys.dopPhaseErr`).** The bracketed number is now `max_i |φ_i(t_app) − 2π f_d,i t_app|`, the difference between the integrated phase and the previous first-order phase, i.e. how much the first-order model would have been off. Typical values (v = 300 km/h, θ₁ = 45°, θ₂ = −30°, τ = 0): d_min = 5 m: 2.4 rad, 30 m: 0.41 rad, 50 m: 0.25 rad, 500 m: 0.025 rad (τ = 10 ms: 13, 2.4, 1.4, 0.15 rad).

## 5. Legacy model (`model = 'legacy'`, default until the tests of `tests/` all pass)

Unchanged: i.i.d. block fading per snapshot, diffuse M = 8 re-drawn each snapshot, 20 % sliding window, expected-covariance SINR, ICI = `K/(K+1)(1−sinc²ε) + 1/(K+1)·⟨1−sinc²(ε_m cosα)⟩_α` (isotropic α). `core.js` is bit-identical to the pre-refactor page (see `tests/legacy_equivalence.js`).

## 6. Training data mode (`trainMode`, both models)

`R̂ = (1/L)Σ x_n x_nᴴ` is built from the snapshot window. Two settings:

| `trainMode` | training snapshots | design name |
|---|---|---|
| `'withSignal'` (default) | target + jammer + noise (the previous behaviour) | **MPDR** (minimum-power distortionless response) |
| `'signalFree'` | jammer + noise only | **MVDR** (textbook assumption: R = R_n) |

Physical basis: with the target in the training data, minimising output power under a distortionless constraint on the *nominal* a(θ̂₁) also minimises the target's own contribution whenever the true target vector h differs from a(θ̂₁) (fading, calibration, aging, angle spread): the filter partly cancels the signal (self-nulling). With signal-free training this mechanism is absent. The signal-free mode is an idealisation (a receiver cannot normally separate the target from its training window); it is the reference for the MPDR loss.

Pairing: in `'signalFree'` the random draws are exactly those of `'withSignal'` (target symbols and diffuse gains are drawn and then not added), so for one seed the interference and noise sequences are identical in both modes. Evaluation (SINR, EVM) always includes the target.

Algorithms: **SMI, DL, BEAMSPACE (and the tapered GSC variants of SMI/DL)** follow `trainMode`. **FOURIER** does not use R̂. **MMSE** always keeps the target: the Wiener filter minimises E|wᴴx − d|² with R = E[xxᴴ] = R_n + P_s h hᴴ and r_xd = P_s a, so R̂ must contain the target for w = R̂⁻¹r_xd to be the Wiener solution; with a signal-free R the same formula gives P_s·R_n⁻¹a, i.e. the MVDR direction, no longer MMSE. In code the target part of each snapshot is stored separately and added back for MMSE. MMSE is therefore unaffected by `trainMode` (and equals SMI-MPDR up to a scalar).

UI: selector "訓練資料" (Rx card, URL flag `?train=signalFree`); names in legends and diagnosis read SMI-MPDR / SMI-MVDR etc. The CSV is unchanged (the `algorithm` field keeps the historical labels SMI-MVDR …). **CSV columns to add in a later round:** `param,train_mode`; `sweep_info` should carry `train=`; (Commit 6/7) `param,gamma_rel_dB`, `param,d_min`, `result,max_angle_drift_deg`.

## 7. Diagonal loading relative to the noise power (DL, unified model only)

Unified: `γ = γ_rel · σ_n²`, UI control `γ_rel` in dB (default **+10 dB**, range −10 … +30 dB), `Sys.gammaRelDb`.
Legacy: unchanged, absolute `γ` (`gammaDL`, default 0.01, log slider 10⁻⁴ … 1); only one of the two controls is shown, according to the model.

**σ_n² definition.** `σ_n² = 10^(−SNR/10)`, the noise power **per array element and per snapshot** (complex, CN(0, σ_n²)). The reference point of the SNR is a *single element's input*: the desired signal has unit power `P_s = 1` at the source and the channel has unit total gain (Σ_i|β_i|² = 1, |Γ_n| = 1), so the per-element signal power is 1 and SNR is the per-element input SNR (array gain is not included). The jammer power is `P_j = 10^(−SIR/10)` in the same units.

Absolute values at the default SNR = 20 dB (σ_n² = 0.01):

| γ_rel | γ |
|---|---|
| 0 dB | 0.01 |
| +10 dB (default) | 0.1 |
| +30 dB | 10 |
| −10 dB | 0.001 |

(For another SNR, γ = 10^((γ_rel − SNR)/10).) Rationale: the loading that matters is relative to the noise floor — the smallest eigenvalues of R̂ are ≈ σ_n² — so γ_rel keeps the same regularisation strength when SNR changes.

**Other constants in the weight designs (listed, not modified):**

| constant | where | absolute or relative |
|---|---|---|
| `LAMBDA_Q = 30` | GSC tapered path: `λ = max(30·σ_n², 10·w_qᴴR̂w_q)` | relative (to σ_n² and to the quiescent output power) |
| `REL_Q = 10` | same | relative |
| `1e-6` pivot clamp in `invertMatrix(..., force)` | SMI/MMSE/BEAMSPACE inversion when a pivot is ~0 | **absolute** (deliberate: it demonstrates the collapse of an exactly singular R̂) |
| `1e-300` test on `a_Bᴴ R_B⁻¹ a_B` and `aᴴR⁻¹a` | weight normalisation guard | absolute (underflow guard only) |
| `1e-12·λ_max` | κ = ∞ decision | relative |
| `1e-3` norm threshold | blocking-matrix Gram–Schmidt | relative (vectors are unit-scale) |
| BEAMSPACE | **no loading constant**: K = 3 DFT beams (nearest to the target + the 2 strongest) | — |

## 8. Carrier frequency, numerology and the ICI diagnosis thresholds

UI (Tx card): `fc` (GHz; presets 3.5 / 5 / 7 / 15 / 28 and a custom field) and the subcarrier spacing Δf (15 / 30 / 60 / 120 kHz, μ = 0…3). They are **independent** (no automatic link between fc and numerology). URL flags `?fc=28&scs=120`. Defaults: fc = 5 GHz, Δf = 15 kHz (the legacy path is bit-identical at these values: T0a).

Everything derived follows `Sys.fc` and `Sys.scs` at every computation: `f_m = v·fc/c`, `ε = f_d/Δf`, `λ = c/fc` (used by the integrated Doppler phase of §4), `T_snap = (1 + cpRatio)/Δf`, the ICI floor (a function of `f_m/Δf`), and the element spacing `d = 0.5λ` (the steering vector depends on `d/λ = 0.5` only, so the physical spacing scales with fc: 30 mm at 5 GHz, 5.4 mm at 28 GHz). `d_min`, v and the angles do not depend on fc.

Consequences worth knowing: ε and the ICI floor depend on `fc/Δf` only; `ρ = L·T_snap·B_D` scales as `fc/Δf`; the aging ratio `τ·B_D` scales with fc and not with Δf (T9).

**ICI diagnosis thresholds** (replace the old `|f_d| > 100 Hz`, which ignored Δf): `CONFIG.iciWarnDb = −30`, `CONFIG.iciSevereDb = −20` on `N_ICI/S` (dB): above −30 dB = warning, above −20 dB = critical. Basis: when ICI is the only impairment `EVM = √(N_ICI/S)`; −30 dB ⇒ EVM 3.2 %, −20 dB ⇒ EVM 10 %. Against the EVM that the existing SER formula needs (SER 10⁻³ / 10⁻⁶): QPSK 30 % / 20 %, 16-QAM 13.1 % / 9.0 %, 64-QAM 6.3 % / 4.4 %. So −30 dB (3.2 %) is the point where ICI starts to eat the margin of the highest-order format (64-QAM needs −27.2 dB for SER 10⁻⁶, and −30 dB is only 2.8 dB inside that), and −20 dB (10 %) is where ICI alone makes 16-QAM miss SER 10⁻⁶ (needs −20.9 dB) and 64-QAM miss 10⁻³ (needs −24.0 dB; SER at −20 dB ≈ 5×10⁻²). The two numbers are round-number choices consistent with those links, not derived optima (**待確認** if a different target SER is wanted).

## 9. Read-outs added during the investigation

- `SINR_opt` (`Sys.sinrOptDb`, panel D; read-only, not used by any algorithm): genie bound `w_opt = R_in⁻¹h`, `SINR_opt = hᴴR_in⁻¹h`, `R_in = P_j g gᴴ + σ²I` (the true interference + noise covariance of the unified model, `g = Γ⊙a(θ₂)`), `h` = this realisation's true channel vector at t_app. Closed form (rank-one inverse): `SINR_opt = (‖h‖² − P_j|gᴴh|²/(σ² + P_j‖g‖²))/σ²`. Legacy model: no per-trial h exists, so the bound is the expected-covariance version, `λ_max(R_in⁻¹ R_t)` (power iteration), with the legacy Rician jammer in `R_in`.
- MMSE's cross-correlation vector: `r_xd = P_s·a(θ̂₁)` — the **nominal** steering vector at the estimated angle (core.js, MMSE branch of `computeMath`); a "model-based" MMSE. A data-driven `r_xd = (1/L)Σ x_n s_nᴴ` was only evaluated in `tests/diag_c10b_mmse.js`; the snapshots now carry the known target symbol (`s1r`, `s1i`) for that purpose. The existing MMSE is unchanged.

## 10. Seeds and CSV schema v2

- UI: seed field + "重新抽種子" button; `?seed=<uint32>`. The seed fixes the phase mismatch (`Core.deriveSeed(seed, 0xCA1, cal_roll_index)`; the "Re-roll φₙ" button increments `cal_roll_index`) and the whole sweep. The live view itself is not frame-reproducible (it depends on timing).
- Sweep: trial k of algorithm a at point i (speed index) uses `Core.deriveSeed(Core.deriveSeed(seed, i), a, k)`; the inner value is the CSV's `point_seed`. `Core.withSeed` restores the live generator afterwards, so a single point can be re-run alone with identical numbers (T10).
- CSV `csv_schema_version = 2` (row right under the header): new `param` rows model, train_mode, gamma_rel_dB, d_min, M, K, sigma_theta, sigma_phi, P_s, tau, theta1, theta2, fc_GHz, delta_f_kHz, seed, cal_roll_index (extra), point_seed (live state: n/a; per point in the sweep table), `result` rows max_angle_drift_deg and SINR_opt_dB, `algorithm` with the new names plus `algorithm_legacy`, `sweep_info` with `train=`, `model=`, `seed=`, and a `point_seed` column in the sweep table. The timestamp is its own row (`meta,timestamp`), the only row that differs between two runs with equal seed and parameters.

## 11. The two Wiener (MMSE) variants

| name (UI, CSV `algorithm`) | key | `r_xd` | CSV `algorithm_legacy` |
|---|---|---|---|
| **MMSE-M** (model-based) | `'MMSE'` | `P_s·a(θ̂₁)` — the nominal steering vector at the estimated angle (unchanged; numbers identical to the former "MMSE", T0a) | `MMSE` |
| **MMSE-P** (pilot-trained) | `'MMSEP'` | `r̂_xd = (1/L)Σ x_n·conj(s_n)`, s_n = the known training symbols (needs known pilots) | `n/a` |

Both solve `w = R̂⁻¹ r_xd` with R̂ containing the target (the Wiener definition: R = R_n + P_s h hᴴ), independent of `trainMode`, and both invert without a safeguard (L < N collapses, see §12). MMSE-P is the only algorithm that uses the known symbols.

Scaling: the SINR does not depend on the scale of w. For EVM, perfect scaling (`g = wᴴh(t_app)`) is what `Sys.evm` uses for every algorithm, MMSE-P included. For MMSE-P (unified) the read-out `Sys.evmPilot` shows the EVM when the output is normalised with the pilot-based channel estimate instead: `ĥ = r̂_xd/P_s`, `ĝ = wᴴĥ`, `ρ = g/ĝ`, `EVM² = |ρ−1|² + |ρ|²(1/SINR + N_ICI/S)`. Effect of the estimation error of ĥ: ĥ is the window average of `h(t_n)`, so (i) the noise/interference cross terms add a random error of order `1/√L`, and (ii) when the channel rotates inside the window (Doppler), `|ĥ|` shrinks (an average of a rotating phasor) and the phase is that of the window centre, not of `t_app`: ρ is then far from 1 and the pilot-scaled EVM is large even if the SINR is high. The same rotation also degrades the direction of `r̂_xd` itself (v = 300 km/h, unified: MMSE-P SINR ≈ −12 dB at K = 20 dB, SINR_opt 29 dB). (**待確認**: the specification asks to normalise EVM with `wᴴĥ`; I kept the headline EVM with perfect scaling and show the pilot-scaled value as a read-out, since it is a second definition.)

## 12. Small snapshots (L < N): what inverts what, and `smiSingular`

Before Commit 15 (the numbers of `tests/diag_c15_paths.js`, N = 8, unified, K = 20 dB, SNR 20 dB, 300 trials):

| algorithm | inversion path | numerical safeguard |
|---|---|---|
| SMI | `mvdrWeights(R̂, a, force = true)` → `invertMatrix(R̂, true)` (Gauss–Jordan, partial pivoting) | a pivot with \|pivot\|² < 1e-12 (**absolute**: \|pivot\| < 1e-6) is replaced by 1e-6 and the elimination goes on: the weights are huge and meaningless ("collapse") |
| DL | same routine with `force = false` on R̂ + γI | none needed: R̂ + γI is non-singular (eigenvalues ≥ γ); if a pivot still failed, the quiescent weight would be used |
| BEAMSPACE | K×K (K = 3) `invertMatrix(R_B, true)` | same absolute clamp, but R_B is singular only for L < K = 3 |
| MMSE-M, MMSE-P | `wienerWeights(R̂, r_xd, true)` → `invertMatrix(R̂, true)` | same absolute clamp as SMI |

R̂ = (1/L)Σ x_n x_nᴴ has rank min(L, N) exactly; the zero eigenvalues are rounding noise (≤ 2·10⁻¹⁴ against λ_max ≈ 80).

| L | rank | κ reported (relative test λ_min ≤ 1e-12 λ_max) | smallest non-zero eigenvalue [min, median, max] | κ on the non-zero part (median) |
|---|---|---|---|---|
| 2 | 2 | ∞ | [1.4e-2, 3.8, 10] | 2.3e1 |
| 4 | 4 | ∞ | [2.1e-3, 8.1e-3, 2.2e-2] | 1.0e4 |
| 6 | 6 | ∞ | [6.6e-5, 1.6e-3, 6.2e-3] | 5.2e4 |
| 8 | 8 | finite | [2.3e-7, 1.5e-4, 1.3e-3] | 5.4e5 |

The "collapse" for L < N is thus produced by the absolute 1e-6 pivot clamp (an implementation artefact, not a physical quantity). At L = N the same absolute threshold can act on a full-rank matrix: λ_min < 1e-6 in 0.2 % of 2000 trials (smallest 3.8e-8).

**`smiSingular`** (SMI only; `CONFIG.smiSingular`, default `'pinv'`): `'clamp'` = the behaviour above (legacy); `'pinv'` = Moore–Penrose pseudo-inverse of R̂ keeping the eigenvalues above `epsRank·λ_max` (`CONFIG.epsRank = 1e-10`, relative), used when the numerical rank is < N: `w = R̂⁺a / (aᴴR̂⁺a)`. The distortionless constraint wᴴa = 1 holds directly by this normalisation (when aᴴR̂⁺a ≠ 0, i.e. a is not orthogonal to the range of R̂; otherwise the quiescent weight is used). Caveat: if a has a component in the null space of R̂ the exact constrained optimum has zero output power from the training data; the pseudo-inverse solution is the one confined to the range of R̂ (the minimum-norm-type choice), not a physically optimal one. With full rank the direct inverse is used (same numbers as before; `pinvHermitian` agrees with it to 1e-12, T12a). UI: status `RANK-DEFICIENT (pinv)` and the message "秩虧：R̂ 秩 = L，求逆為偽逆" replace the "collapse" message of SMI. The legacy-equivalence test T0a runs with `smiSingular = 'clamp'`. DL, BEAMSPACE and MMSE are unchanged (DL does not read the option, T12d). Condition number read-out: `Sys.kappaRank` = λ_max/λ_min, shown as ∞ with the rank (`rank/N`) when the numerical rank (`Sys.rankR`, threshold `epsRank·λ_max`) is below N.

## 13. Modulation in the signal chain and symbol-level EVM

**Target symbols follow the modulation.** The target symbol s_n of the training snapshots (both models) is drawn from the selected constellation (QPSK / 16-QAM / 64-QAM, unit average power; `Sys.txSym`); with QPSK it uses the original draws, so the results are bit-identical to before (T0a runs with `mod = 'QPSK'`, the only symbol set the old page used). The UI default modulation is 16-QAM, so the default numbers of the page change slightly (non-constant modulus: |s_n|² fluctuates, the sample covariance has a different finite-sample noise; the expected covariance is the same). The jammer symbol stays QPSK. Tests that need the numbers of the previous round (T11a) select QPSK explicitly.

**Symbol-level EVM (unified only; `Core.symbolLevel`, `Sys.symbolEvm(Ns = 4000, keep)`).** For one channel realisation with output gain g = wᴴh(t_app) (`Sys.gR/gI`, S = |g|²): a fresh random symbol s (selected modulation) is sent Ns times; the receiver output after gain normalisation, assuming **a perfect channel estimate (g known exactly)**, is

```
ŝ = s + ( wᴴj + wᴴn + ICI ) / g
```

with wᴴn ~ CN(0, Nn), wᴴj = √I × (unit-power QPSK interferer symbol, random phase), ICI ~ CN(0, S·N_ICI/S): **a Gaussian approximation of the OFDM leakage** (the true ICI of a few strong paths is not Gaussian). EVM_meas = √(Σ|ŝ − s|² / Σ|s|²). The first 500 outputs are the constellation scatter of the unified model (the legacy scatter remains synthetic: EVM-scaled Gaussian). Panel C shows the measured EVM next to the analytic one. This does not change `Sys.evm` (analytic, used by the sweep and the SER).

**Sweep.** Option "符元層級驗證" (default off): per sweep point and algorithm, 5 further channel realisations × 4000 symbols (child seeds 1000 + j); mean EVM → CSV column `evm_meas_sample_<alg>_pct` (n/a when off; `sweep_info` carries `symbol_validation` and `evm_meas_sample_n`).

Limits (T13 header): the measured chain is built from the same S, I, N_n, N_ICI/S as the analytic EVM, so T13a verifies the arithmetic of the chain (scaling, normalisation, estimators), not the ICI physics.
