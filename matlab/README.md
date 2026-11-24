# MATLAB / Simulink sources

Files from our group's SSY285 submissions, copied unchanged. They were **not run** when building this repo (no MATLAB
available); `../analysis/` re-implements the same model in Python and checks it against the numbers in the reports.

| Folder | Files | What it does |
|---|---|---|
| `E1/` | `E1C.m` | state-space model, Laplace-domain steady state, `lsim` check |
| | `E1D.m`, `E1E.m` | symbolic Jacobian linearisation, step and 5 Hz responses, five initial load voltages |
| | `E1F_new.m` | transfer function, `zpkdata`, `minreal` |
| `E2/` | `assignment_E2.mlx`, `a_1.m` | symbolic controllability / observability, PBH tests, discretisation |
| `E3/` | `E3.mlx` | Kalman filter, `dlqr`, `lqgreg` script |
| | `E3.slx` | Simulink LQG loop with process and measurement noise |

Provenance: the course template files and lecture material are not included. The files above are the group's own work product;
the `.slx` and `.mlx` metadata carry the OS username of the account that saved them.
`E3.mlx` uses `Qx = 50 I, Qu = 0.01 I` while the E3 report prints results for `Qx = diag(0,0,0,1e4), Qu = 0.001`
(see `../docs/results.md`).
