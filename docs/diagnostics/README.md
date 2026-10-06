# Diagnostics (informational; not part of `node tests/run_all.js`)

Produced by `tests/diag_b1_mmsep_L.js`, `tests/diag_b2_t13a.js` (`--only-blocks`), `tests/diag_b4_afc.js`; B3 is the informational block T17b2 in `tests/t17_pointing.js`.
No criterion, tolerance or algorithm was changed for any of them.

| file | question | observation |
|---|---|---|
| `b1_mmsep_L.txt` | Is the MMSE-P jump at L = 8 -> 9 (T11c / T16c) tied to an ill-conditioned R_hat at L ~ N? | The mean SINR of MMSE-P (17.44 dB) and of SMI signalFree (21.43 dB) has its minimum exactly at L = N = 8, where the median condition number of R_hat is 5.5e5 (1.9e5 at L = 9, 1.2e5 at L = 10), the median smallest non-zero eigenvalue is smallest (1.4e-4) and the median weight norm largest (1.12 against 0.75 at L = 9). The "jump" is the recovery from this L = N valley. Direct inverse and forced pseudo-inverse give the same SINR (difference <= 1e-5 dB), so the inversion routine is not the cause. |
| `b2_t13a_ablation.txt`, `b2_t13a_blocks.txt` | Why is T13a (signalFree, v = 300) at +0.20 % (3.72 SE)? | The change from +0.07 % to +0.20 % appears at Commit 19 (the jammer waveform default changed the random numbers; every earlier commit gives +0.0708 %, every later one +0.1957 %). In ten consecutive blocks of 300 realisations only block 0 (the test seeds) is high (z = 3.72); the other nine have z between -1.45 and 1.10, the mean of the ten is +0.022 % +- 0.023 %. With 2000 realisations: +0.033 % +- 0.020 %. With the ICI term switched off the same seeds give -0.019 %. No bias is found; the failing value is a fluctuation of that seed set. |
| T17b2 (in the T17 output) | Is the finite-L deviation of T17b the Reed-Mallett-Brennan loss? | see the T17 summary table |
| `b4_afc.txt` | Does a periodogram frequency estimate repair MMSE-P? | See the file; with the criteria of the earlier T15 specification the periodogram variant fails T15a (v = 0: -1.93 dB, |z| = 6.2) and, at theta1 = 45, T15b (genie - on = 4.05 dB). Not implemented. |

## Sensitivity of E6 (S1-S3; `tests/diag_s_sensitivity.js s1|s2|s3`, 500 realisations per cell)

| file | observation |
|---|---|
| `s1_sensitivity.csv/.txt` | 240 cells ((fc, df) x d_min x tau x v x {SMI, DL}, signalFree): the aging loss reaches the ICI loss in 0 cells. Closest: 28 GHz / 120 kHz, d_min 10 m, tau 10 ms, 300 km/h, SMI: aging 5.53 dB, ICI 9.91 dB (ratio 0.56). E6 itself uses 5 GHz / 15 kHz (aging 1.20 dB, ICI 12.68 dB for SMI at 300 km/h). |
| `s2_dl_vs_smi.txt` | At 300 km/h DL loses 4.91 dB, SMI 1.35 dB (gamma_rel +10 dB). The DL loss falls to the SMI value for gamma_rel <= -4 dB (best at -4 dB: 27.35 dB, loss 1.39 dB) and grows with gamma_rel (7.7 dB at +24 dB); at v = 0 the loading gains 0.27 dB at +10 dB. Gain at the jammer direction: at the estimation time -46.3 dB (DL) against -54.9 dB (SMI), at the application time -35.7 dB against -44.8 dB (jammer drift 1.15 deg, +10.1 / +10.6 dB). Null width (below -30 dB) is about the same (11.1 deg DL, 10.0 deg SMI at v = 300). So the extra DL loss goes with the shallower null (8.7 dB at the estimation time, 9.1 dB at the application time), not with a narrower one. |
| `s3_withsignal.txt` | SMI withSignal: with K = inf (no diffuse paths) the SINR does not rise with speed (11.4 dB at v = 0, 6.0 dB at 200 km/h, 6.7 dB at 300 km/h with tau = 10 ms; 6.1 dB with tau = 0); with the original K = 20 dB it rises from -5.5 dB to 7.2 dB (125 km/h) and falls to 4.6 dB; diffuse only: from -21.8 dB to 16.6 dB (75 km/h) and falls to 8.7 dB. The rise is therefore due to the diffuse paths (Doppler differences between them). The share of the first three eigenvalues of the target covariance stays at 99.7-100 % for K = 20 dB (diffuse only: 79.2 % / 96.9 % / 99.8 % for the 1st / 1st+2nd / 1st+2nd+3rd eigenvalue at 225 km/h and 75.6 % / 95.4 % / 99.6 % at 300 km/h, from 100 % at v = 0), i.e. the rise is not visible as a rank change of the target covariance at K = 20 dB. |

## L < N: MMSE-P in `verify_cases.m` (case `static_L4_signalFree_d3`)

Scripts (read-only on the repository; run from anywhere with `node docs/diagnostics/<file>.js`, outputs stored next to them as `.txt`):

| script | what it does | result |
|---|---|---|
| `l4_mmsep_vs_qr.js` | for the four L = 4 cases of `data/matlab_cases.json`: eigenvalues of the Wiener R_hat (kept / dropped by the `epsRank` threshold), kappa of the retained part, and the web weights against a QR (Gram-Schmidt, applied twice) reference of the same minimum-norm solution `w = (X^H)^+ conj(s)` | retained kappa 4.1e3 - 1.9e4 (the full-matrix kappa of about 1e17 is rounding noise of the zero eigenvalues); smallest kept eigenvalue 5e-5 ... 2.5e-4 of lambda_max, largest dropped 1e-16: no eigenvalue is near the threshold 1e-10. Web weights against the QR reference: relative error 2.8e-10 / 5.7e-10 / 6.6e-10 / 4.8e-10 (SINR 3.6e-10 / 2.1e-9 / 6.4e-10 / 2.2e-10 dB). For `static_L4_signalFree_d3` these are the numbers (5.72e-10 and 2.13e-9 dB) that `verify_cases.m` reported; MATLAB itself was not run here |
| `l4_mmsep_jacobi_and_residual.js` | Jacobi stopping tolerance of `pinvHermitian`, and the residual of the consistent system `X^H w = conj(s)` (exact solution: 0) | tolerance contributes 3.4e-12 only; residual 1.9e-11 for the pseudo-inverse solution against 2.4e-16 for the QR solution; reordering the sum of R_hat (a 1.2e-16 relative change) changes the pseudo-inverse weights by 7.4e-10 and the SINR by 1.5e-9 dB |
| `l4_mmsep_rounding_sensitivity.js` | 400 random L = 4 realisations (unified, signalFree, K = 20 dB): relative change of the pseudo-inverse weights when the summation order of R_hat is changed, against kappa of the retained part | median 2.5e-10, 90th percentile 1.4e-9, maximum 1.2e-8; 14 % of the realisations exceed 1e-9; log-log slope 1.7 against kappa (correlation 0.81) |

Reading: for L < N the pseudo-inverse solve is sensitive to rounding of R_hat; the tolerance 1e-9 of `verify_cases.m` is tight for these cases. The tolerance was not changed; `verify_cases.m` lists the rank-deficient cases separately.

