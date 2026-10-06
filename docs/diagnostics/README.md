# Diagnostics (informational; not part of `node tests/run_all.js`)

Produced by `tests/diag_b1_mmsep_L.js`, `tests/diag_b2_t13a.js` (`--only-blocks`), `tests/diag_b4_afc.js`; B3 is the informational block T17b2 in `tests/t17_pointing.js`.
No criterion, tolerance or algorithm was changed for any of them.

| file | question | observation |
|---|---|---|
| `b1_mmsep_L.txt` | Is the MMSE-P jump at L = 8 -> 9 (T11c / T16c) tied to an ill-conditioned R_hat at L ~ N? | The mean SINR of MMSE-P (17.44 dB) and of SMI signalFree (21.43 dB) has its minimum exactly at L = N = 8, where the median condition number of R_hat is 5.5e5 (1.9e5 at L = 9, 1.2e5 at L = 10), the median smallest non-zero eigenvalue is smallest (1.4e-4) and the median weight norm largest (1.12 against 0.75 at L = 9). The "jump" is the recovery from this L = N valley. Direct inverse and forced pseudo-inverse give the same SINR (difference <= 1e-5 dB), so the inversion routine is not the cause. |
| `b2_t13a_ablation.txt`, `b2_t13a_blocks.txt` | Why is T13a (signalFree, v = 300) at +0.20 % (3.72 SE)? | The change from +0.07 % to +0.20 % appears at Commit 19 (the jammer waveform default changed the random numbers; every earlier commit gives +0.0708 %, every later one +0.1957 %). In ten consecutive blocks of 300 realisations only block 0 (the test seeds) is high (z = 3.72); the other nine have z between -1.45 and 1.10, the mean of the ten is +0.022 % +- 0.023 %. With 2000 realisations: +0.033 % +- 0.020 %. With the ICI term switched off the same seeds give -0.019 %. No bias is found; the failing value is a fluctuation of that seed set. |
| T17b2 (in the T17 output) | Is the finite-L deviation of T17b the Reed-Mallett-Brennan loss? | see the T17 summary table |
| `b4_afc.txt` | Does a periodogram frequency estimate repair MMSE-P? | See the file; with the criteria of the earlier T15 specification the periodogram variant fails T15a (v = 0: -1.93 dB, |z| = 6.2) and, at theta1 = 45, T15b (genie - on = 4.05 dB). Not implemented. |
