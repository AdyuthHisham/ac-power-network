# Results: data sources and mismatches

All numbers in the README come from the scripts in `analysis/`, which write JSON to `docs/metrics/`.
`docs/report_comparison.md` (generated) lists every number printed in the three SSY285 reports next to the value computed here.

## Data sources

| README item | Script | Output | Source of the number |
|---|---|---|---|
| Circuit diagram | `analysis/diagram.py` | `media/network.svg` | parameters from the E1 report and the `.m` files |
| Steady-state amplitude and phase | `analysis/e1_model.py` | `media/e1_steady_state.png`, `metrics/e1.json` | phasor solve of the 4-state model, checked against a time simulation |
| Step response, initial-condition response | `analysis/e1_model.py` | `media/e1_step.png`, `media/e1_initial.png` | `scipy.signal.lsim` |
| Poles, zeros, transfer functions | `analysis/e1_model.py` | `media/e1_pzmap.png`, `metrics/e1.json` | `ss2tf`, `eig` |
| Unobservable load values | `analysis/e2_structure.py` | `media/e2_pbh.png`, `metrics/e2.json` | sympy: `det(Obs) = k (150 LL - RL)(300 LL - RL)/LL^2` |
| Condition numbers, discrete eigenvalues | `analysis/e2_structure.py` | `media/e2_discrete.png`, `metrics/e2.json` | matrix exponential (ZOH, Ts = 5 ms) |
| Kalman gain, LQR gain, estimation and tracking numbers | `analysis/e3_lqg.py` | `media/e3_*.png`, `media/hero_lqg.gif`, `metrics/e3.json` | discrete Riccati equations; 100 Monte Carlo runs, seeds 1000-1099 |
| Report-vs-repo table | `analysis/report_check.py` | `docs/report_comparison.md` | report values transcribed in that script |
| Simulink diagram | cropped from the E3 report, page 6 | `media/e3_simulink_from_report.png` | `pdftoppm -r 200 -f 6 -l 6 -x 240 -y 540 -W 1200 -H 800 -png -singlefile SSY285_ProjectE3.pdf e3_simulink_from_report` |

Everything except the Simulink crop is a **Python re-implementation** of the MATLAB model. Nothing was run in MATLAB or Simulink here.

## Where the Python numbers and the reports differ

1. **`Bd`, fourth row.** Report: 0.4932 / 0.3824. ZOH (`expm`) gives 0.4956 / 0.3837. Rows 1-3 and all of `Ad` agree to four decimals. Effect: discrete controllability condition number 439.25 here vs 437.857 in the report. Cause not identified.
2. **LQR gain.** With the report's weights (`Qu = 0.001`, `Qx = diag(0,0,0,1e4)`) the gain is 1.5 % or less away from the printed one (e.g. 120.20 vs 119.79). Cause not identified; using the 4-decimal `Bd` from the report does not reproduce it.
3. **Kalman gain naming.** The report's `K` equals the innovation gain `M = P Cd'(Cd P Cd' + R2)^-1` (matches to under 0.001). The observer eigenvalues printed in the report are those of `Ad - M Cd`. The error dynamics of the predictor form are `Ad - Ad M Cd`: -0.0425 +/- 0.4289i, 0.3927, 0.1571. Both sets lie inside the unit circle, so the stability statement holds either way. The feedback loop in the Python simulation uses the predictor form.
4. **`P`.** Largest entry difference is 0.7 % of max |P| (`P(4,4)`: 9.3774e4 report, 9.4428e4 here).
5. **`R2(1,1)`.** The report prints 4878.78; its own variance line and `E3.mlx` give 48787.78, which is what is used.
6. **`E3.mlx` vs the report.** `E3.mlx` designs the LQR with `Qx = 50 I, Qu = 0.01 I`; the report prints results for the other weights. Both were simulated: the noise-free 5 Hz tracking is the same to the digits in the table below. The MATLAB output text stored in `E3.mlx` (observer eigenvalues -0.0419 +/- 0.4269i, 0.3888, 0.1619) is reproduced to four decimals when the 4-decimal `Ad`, `Bd` typed in that file are used.
7. **Feed-forward gain `K_FF`.** The report's Simulink diagram has a `K_FF` block but gives no value. The Python simulation uses the minimum-norm DC feed-forward `pinv(Cx (I - Ad + Bd L)^-1 Bd)`, value in `metrics/e3.json`.
8. **Linearisation.** The plant is already linear, so the Jacobian linearisation at zero returns `A`, `B`, `C` unchanged (finite-difference check, `metrics/e1.json`).
9. **Pole listing.** The report prints the poles with a common factor 1e3 (-0.1202 +/- 2.8248i, -0.2035, -0.1061). That reads as -120.2 +/- 2824.8i, -203.5, -106.1 and matches.

## Estimation and tracking results (Python re-implementation)

Setup: report parameters, `R1 = 222222.2 I`, `R2` as in the report, Ts = 5 ms, reference `1000 cos(2 pi 5 t)` V, 10 s runs, errors measured for t > 1 s.

| Quantity | Value |
|---|---|
| Raw load-voltage sensor error (RMS, 100 runs) | 220.5 V |
| Kalman estimate, predicted `x[k|k-1]` (the state the controller uses), RMS error | 306.9 V (theory 307.3 V) |
| Kalman estimate, measurement-updated `x[k|k]`, RMS error | 153.8 V (theory 153.6 V) |
| Tracking error with noise, RMS | 335.2 V |
| Noise-free 5 Hz tracking | gain 1.000, lag 9.0 deg (= one 5 ms sample), RMS error 111.0 V |

The predicted estimate is worse than the sensor, as the theory says it should be, because the process-noise model (`Nd = Bd`,
sigma 471 V on each input) is large. The updated estimate is 30 % better than the sensor.

## Not done

- No tracking-error metric exists in the reports (the Simulink scopes are qualitative). The numbers above are from the Python simulation only.
- No hardware, no real measurements. The inverted-pendulum lab in the source folder is course material and is not part of this repo.
- Simulink and MATLAB were not run.
