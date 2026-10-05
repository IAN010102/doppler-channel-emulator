# PARAMS.md — parameters, conventions and the two channel models

`core.js` holds the DOM-free physics (browser global `Core`, or `require('./core.js')` in Node); `index.html` is UI + drawing only.
All defaults live in `CONFIG` (top of `core.js`); all randomness goes through `Core.rng()` (mulberry32, `Core.setSeed(s)`).
The UI picks a fresh random seed at every page load.

## 1. CONFIG defaults

| key | value | meaning |
|---|---|---|
| `fc` | 5 GHz | carrier frequency |
| `c` | 3×10⁸ m/s | speed of light |
| `scs` | 15 kHz | OFDM subcarrier spacing Δf |
| `N` | 8 | default number of array elements (UI: 2–16) |
| `d_lambda` | 0.5 | element spacing in wavelengths (ULA) |
| `cpRatio` | 0.07 | CP ratio; snapshot period `T_snap = (1 + cpRatio)/Δf` = 71.33 µs (5G-NR-like symbol length) |
| `d_min` | 30 m | perpendicular distance from the ground point to the straight track (§4); renamed from `R_min` |
| `SIGMA_ANG_DEG` | 10° | angular spread σ_θ of the diffuse paths |
| `M_SCAT` | 8 | diffuse paths per snapshot — **legacy model only** |
| `M_UNIFIED` | 32 | diffuse paths per trial — **unified model** |
| `LAMBDA_Q`, `REL_Q` | 30, 10 | loading rule of the tapered (GSC) adaptive path (unchanged) |
| `REFRESH` | 0.2 | sliding-window refresh fraction — **legacy model only** |
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
Reference angles are the slider values = angles **at t_app**; the path angle follows the exact straight-track law θ_i(t) = trackAngle(θ_i,0, v, d_min, t − t_app), θ̇ = v sinθ|sinθ|/d_min for target paths and jammer alike (§4). The Doppler frequency `f_d,i = f_m cosθ_i,0` (angle at t_app) is held constant over the trial: a first-order approximation whose validity is shown by the read-outs of §4.

Channel and snapshots:

```
h(t)  = Σ_i β_i · Γ⊙a(θ_i(t)) · exp(j2π f_d,i t)
x_n   = h(t_n)·s_n + √P_j · Γ⊙a(θ₂(t_n))·exp(j2π f_d,2 t_n)·j_n + noise_n        (s_n, j_n unit-power QPSK, noise CN(0,σ²))
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

**First-order Doppler (kept, with its validity condition).** `f_d,i = f_m cosθ_i,0` is evaluated at the angle at t_app and held constant while the array response follows θ_i(t). The true Doppler phase is ∫f_d dt, so the constant-f_d model makes a phase error `≈ π·|ḟ_d|·T²` on a path (T = t_app = window + τ, the quadratic term relative to t_app). The approximation is valid when this is ≪ 1 rad, i.e. when the angle change within the window (and τ) is small. Diagnostics in panel D ("窗內角度漂移 Δθ_max", the worst phase error in brackets): `Sys.angDrift` = largest |θ_i(first snapshot) − θ_i(last snapshot)| over the target paths and the jammer; `Sys.dopPhaseErr` = the worst π|ḟ_d|T².

Measured (v = 300 km/h, θ₁ = 45°, θ₂ = −30°, L = 100, N irrelevant): 

| d_min | Δθ_max (τ = 0) | phase error τ = 0 | phase error τ = 10 ms |
|---|---|---|---|
| 5 m | 4.99° | 2.5 rad | 14 rad |
| 30 m (default) | 0.87° | 0.42 rad | 2.4 rad |
| 50 m | 0.52° | 0.25 rad | 1.5 rad |
| 500 m | 0.05° | 0.025 rad | 0.15 rad |

So at the default d_min = 30 m the first-order Doppler is only marginally valid at τ = 0 and not valid for τ of several ms (**待確認**: whether to evaluate f_d,i along the trajectory, i.e. use the phase 2π∫f_d dt, instead; not done — outside the agreed scope).

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
