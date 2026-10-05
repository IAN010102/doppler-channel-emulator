# Golden snapshot change log

## Commit 19 — jamWave default 'gaussian' (the jammer symbols draw different random numbers; FOURIER does not depend on the jammer waveform in its weights, so only its target term... unchanged)

Mean-free per-seed SINR [dB], seeds 11/22/33/44, old -> new (the EVM, SINR_opt and SER values in the JSON change accordingly):

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
