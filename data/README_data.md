# data/ — files for the MATLAB scripts and the report figures

All files are produced by Node scripts in `tests/` (no random numbers are drawn when the files are *read*). Regenerate with
`node tests/make_matlab_cases.js` and `node tests/make_curves.js`. Every CSV starts with `#` comment lines (schema version, generator, seed, trials, parameters) followed by the column header.
MATLAB: `readtable(file, 'CommentStyle', '#')`.

Conventions: angles in degrees, SINR in dB, EVM in percent. "mean" is the mean over channel realisations (trials) and "se" the standard error of that mean. Trial *t* of every point uses the seed
*seed0 + t* (the same seeds for all points and algorithms, i.e. the points are paired). Unified model, N = 8, K = 20 dB, SNR = 20 dB (per element), SIR = −10 dB, 16-QAM unless a file says otherwise.

## matlab_cases.json

22 fixed cases (schema_version 1): 20 static cases (v = 0; L ∈ {4, 8, 12, 24, 100}; trainMode signalFree / withSignal; δθ ∈ {0, 3}°) and 2 dynamic cases (v = 300 km/h, τ = 2 ms, d_min = 30 m, L = 100,
both training modes). Complex numbers are `{re, im}` pairs of real arrays; matrices are arrays of rows. The header fields `definitions`, `weight_normalisation` and `pseudo_inverse` give the exact
formulas (R̂, SINR, the normalisation of each algorithm's weights). Per case: `params`, `X` (N×L snapshots, the data that contain the target), `X_sf` (the same draws without the target term),
`s` (known target symbols), `h` (true channel at the application time), `a_J` (jammer steering vector), `a_assumed` (a(θ̂₁+δθ)), `gamma` (hardware phase mismatch), `R_hat_train_web`,
`beamspace` (`bins`, `B`) and `algorithms` (for FOURIER, MMSE-M, MMSE-P, SMI, DL, BEAMSPACE: `w`, `sinr_dB`, `sinr_opt_dB`, `status`, `rank_R`).

## curve_evm_vs_velocity.csv

| column | meaning |
|---|---|
| `fc_GHz` | carrier frequency (5 or 28) |
| `delta_f_kHz` | subcarrier spacing (15 with 5 GHz, 120 with 28 GHz) |
| `velocity_kmh` | speed, 0 … 300 in steps of 10 |
| `algorithm` | FOURIER, MMSE-M, MMSE-P, SMI, DL, BEAMSPACE |
| `evm_mean_pct` | mean analytic EVM = √(1/SINR + N_ICI/S) over the trials, percent |
| `evm_se_pct` | standard error of that mean |
| `ici_floor_evm_pct` | EVM floor caused by ICI alone (no spatial filtering), percent |
| `ici_floor_se_pct` | its standard error (0: it does not depend on the realisation of the noise) |
| `n_trials` | number of channel realisations per point (500) |

## curve_sinr_vs_L.csv

| column | meaning |
|---|---|
| `train_mode` | signalFree (training data without the target) or withSignal |
| `L` | number of snapshots, 2 … 48 |
| `algorithm` | SMI, DL, BEAMSPACE, MMSE-P, or SINR_opt (genie upper bound) |
| `sinr_mean_dB` | mean SINR over the trials |
| `sinr_se_dB` | standard error of that mean |
| `n_trials` | number of channel realisations per point (1000) |

## curve_sinr_vs_delta_theta.csv

| column | meaning |
|---|---|
| `train_mode` | signalFree or withSignal |
| `L` | number of snapshots (100 or 12) |
| `delta_theta_deg` | pointing error δθ, −5 … +5 in steps of 0.5 |
| `algorithm` | SMI, DL, BEAMSPACE, MMSE-M, MMSE-P (not affected by δθ), or SINR_opt |
| `sinr_mean_dB` | mean SINR over the trials |
| `sinr_se_dB` | standard error of that mean |
| `n_trials` | number of channel realisations per point (1000) |

## curve_doa_spectrum.csv

| column | meaning |
|---|---|
| `angle_deg` | angle of the spectrum grid, −90 … 90 in steps of 0.5 |
| `capon_single_dB` | Capon spectrum of trial 0, dB relative to its maximum |
| `music_single_dB` | MUSIC spectrum (2 assumed sources) of trial 0, dB relative to its maximum |
| `capon_mean_dB` | mean over the trials of the Capon spectrum in dB |
| `capon_se_dB` | standard error of that mean |
| `music_mean_dB` | mean over the trials of the MUSIC spectrum in dB |
| `music_se_dB` | standard error of that mean |
| `n_trials` | number of channel realisations (1000) |

Scenario: withSignal training data, θ₁ = −10°, θ₂ = 40°, L = 100 (the spectrum contains target and jammer).
