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
| `R_min` | 30 m | **closest-approach distance to the trajectory** (see §4) |
| `SIGMA_ANG_DEG` | 10° | angular spread σ_θ of the diffuse paths |
| `M_SCAT` | 8 | diffuse paths per snapshot — **legacy model only** |
| `M_UNIFIED` | 32 | diffuse paths per trial — **unified model** |
| `LAMBDA_Q`, `REL_Q` | 30, 10 | loading rule of the tapered (GSC) adaptive path (unchanged) |
| `REFRESH` | 0.2 | sliding-window refresh fraction — **legacy model only** |
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
Reference angles are the slider values = angles **at t_app**; path angle `θ_i(t) = θ_i,0 + θ̇_i·(t − t_app)`, `θ̇_i = −v·sinθ_i,0/R_min` (jammer: `+v·sinθ₂/R_min`, existing sign). The Doppler frequency `f_d,i = f_m cosθ_i,0` is held constant over the trial (first-order consistent with the drift).

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

## 4. About R_min

`R_min` is documented as the closest-approach distance to the trajectory. The drift formula `θ̇ = −v sinθ / R_min` is kept as it was ("existing implementation and sign"). Strictly, for a straight trajectory with closest-approach distance `d`, `θ̇ = −v sin²θ / d`; the kept formula equals it only at θ = ±90° and overestimates the drift by 1/|sinθ| elsewhere. (Flagged, not changed.)

## 5. Legacy model (`model = 'legacy'`, default until the tests of `tests/` all pass)

Unchanged: i.i.d. block fading per snapshot, diffuse M = 8 re-drawn each snapshot, 20 % sliding window, expected-covariance SINR, ICI = `K/(K+1)(1−sinc²ε) + 1/(K+1)·⟨1−sinc²(ε_m cosα)⟩_α` (isotropic α). `core.js` is bit-identical to the pre-refactor page (see `tests/legacy_equivalence.js`).
