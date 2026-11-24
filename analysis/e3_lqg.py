"""E3: discrete Kalman filter + LQR (= LQG) tracking a 5 Hz, 1 kV reference on the discretised network.

Python re-implementation of the MATLAB/Simulink design. The Simulink feed-forward gain K_FF is not
given in the report; here it is the minimum-norm DC feed-forward (see ff_gain)."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from common import *

style()
A, B, C = plant(); Ad, Bd = zoh(A, B)
Cd = np.array([[0, 0, 0, 1.0], [0, 0, 1, 0]])               # measured: load voltage, load current
R1 = 222222.2 * np.eye(2)                                   # process noise on the inputs (Nd = Bd)
R2 = np.array([[48787.78, 73.99], [73.99, 0.4489]])         # measurement noise (report values)
K, P = kalman_gain(Ad, Bd, Cd, R1, R2)                      # predictor-form gain
M = P @ Cd.T @ np.linalg.inv(Cd @ P @ Cd.T + R2)            # innovation gain (what the report prints as K)
Lq = dlqr(Ad, Bd, np.diag([0, 0, 0, 1e4]), 1e-3 * np.eye(2))
Lq_alt = dlqr(Ad, Bd, 50 * np.eye(4), 1e-2 * np.eye(2))     # weights found in E3.mlx

# --- compare with the numbers printed in the report ---
rep = dict(M=np.array([[0, 1.4317], [0.0015, -0.4666], [0.0005, 0.6495], [0.3082, 115.703]]),
           P=1e4 * np.array([[.0062, -.0057, .0003, .0620], [-.0057, .0056, 0, .0015], [.0003, 0, .0002, .0359], [.0620, .0015, .0359, 9.3774]]),
           L=np.array([[119.7861, 123.1789, -126.5657, 0.0358], [94.4463, 98.8154, -100.8499, 0.0326]]))
cmp_ = dict(
    M_max_abs_diff=np.abs(M - rep["M"]).max(),
    P_max_rel_diff=float(np.max(np.abs(P - rep["P"]) / np.abs(P).max())),
    L_max_rel_diff=float(np.max(np.abs(Lq - rep["L"]) / np.abs(rep["L"]))),
    observer_eig_report_formula=np.linalg.eigvals(Ad - M @ Cd),           # = the eigenvalues printed in the report
    observer_eig_error_dynamics=np.linalg.eigvals(Ad - K @ Cd),
    lqr_eig=np.linalg.eigvals(Ad - Bd @ Lq))


def ff_gain(L):
    """Min-norm feed-forward so that the DC gain from r to u_load is 1 (K_FF is not given in the report)."""
    n = Ad.shape[0]
    g = (np.array([[0, 0, 0, 1.0]]) @ np.linalg.solve(np.eye(n) - Ad + Bd @ L, Bd))
    return np.linalg.pinv(g)


def simulate(L, seed=None, T=10.0, noise=True):
    rng = np.random.default_rng(seed)
    n = int(round(T / TS)) + 1
    V = rng.standard_normal((n, 2)) @ np.linalg.cholesky(R1).T if noise else np.zeros((n, 2))
    Wn = rng.standard_normal((n, 2)) @ np.linalg.cholesky(R2).T if noise else np.zeros((n, 2))
    kff = ff_gain(L)
    t = np.arange(n) * TS
    r = 1000 * np.cos(2 * np.pi * 5 * t)
    x, xh = np.zeros(4), np.zeros(4)
    X, XH, XC, Y, U = np.zeros((n, 4)), np.zeros((n, 4)), np.zeros((n, 4)), np.zeros((n, 2)), np.zeros((n, 2))
    for k in range(n):
        v, w = V[k], Wn[k]
        y = Cd @ x + w
        u = (kff * r[k]).ravel() - L @ xh
        X[k], XH[k], Y[k], U[k] = x, xh, y, u
        XC[k] = xh + M @ (y - Cd @ xh)          # measurement-updated (current) estimate, not fed back
        x = Ad @ x + Bd @ (u + v)
        xh = Ad @ xh + Bd @ u + K @ (y - Cd @ xh)
    return t, r, X, XH, XC, Y, U


def rms(a):
    return float(np.sqrt(np.mean(np.square(a))))


def fit_5hz(t, y, r):
    """Gain and phase lag (deg) of y relative to r at 5 Hz, from a least-squares fit on t > 1 s."""
    s = t > 1.0
    Mx = np.column_stack([np.cos(2 * np.pi * 5 * t[s]), np.sin(2 * np.pi * 5 * t[s])])
    a, b = np.linalg.lstsq(Mx, y[s], rcond=None)[0]
    return float(np.hypot(a, b) / 1000), float(np.degrees(np.arctan2(b, a)))


# ---------- deterministic tracking (no noise) for both weightings ----------
res = {}
for name, L in (("report_weights", Lq), ("mlx_weights", Lq_alt)):
    t, r, X, XH, XC, Y, U = simulate(L, noise=False)
    gain, phase = fit_5hz(t, X[:, 3], r)
    res[name] = dict(noise_free_gain_5Hz=gain, noise_free_phase_deg=phase,
                     noise_free_rms_tracking_error_V=rms((X[:, 3] - r)[t > 1]), peak_u_V=float(np.abs(U).max()))

# ---------- Monte Carlo: estimator performance with the report weights ----------
N_MC = 100
est_err, cur_err, raw_err, trk_err = [], [], [], []
for s in range(N_MC):
    t, r, X, XH, XC, Y, U = simulate(Lq, seed=1000 + s)
    m = t > 1.0
    est_err.append(rms((XH[:, 3] - X[:, 3])[m])); cur_err.append(rms((XC[:, 3] - X[:, 3])[m])); raw_err.append(rms((Y[:, 0] - X[:, 3])[m]))
    trk_err.append(rms((X[:, 3] - r)[m]))
res["monte_carlo"] = dict(n_runs=N_MC, seeds="1000..1099", rms_load_voltage_predicted_estimate_error_V=float(np.mean(est_err)),
                          rms_load_voltage_updated_estimate_error_V=float(np.mean(cur_err)),
                          rms_load_voltage_raw_measurement_error_V=float(np.mean(raw_err)),
                          rms_tracking_error_V=float(np.mean(trk_err)),
                          theory_posterior_std_load_voltage_V=float(np.sqrt(
                              (P - M @ Cd @ P)[3, 3])), theory_prior_std_load_voltage_V=float(np.sqrt(P[3, 3])))
save_metrics("e3", dict(design_vs_report=cmp_, results=res, K_predictor=K, M_innovation=M, P=P, L_report_weights=Lq,
                        L_mlx_weights=Lq_alt, K_FF=ff_gain(Lq)))

# ---------- figures ----------
t, r, X, XH, XC, Y, U = simulate(Lq, seed=7)
m = t <= 2.0
fig, ax = plt.subplots(2, 1, figsize=(7.5, 5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.4]))
ax[0].plot(t[m], Y[m, 0], color=OI["sky"], lw=0.8, label="noisy measurement")
ax[0].plot(t[m], X[m, 3], color=OI["blue"], lw=1.4, label="true load voltage")
ax[0].plot(t[m], XC[m, 3], color=OI["orange"], lw=1.4, label="Kalman estimate (updated)")
ax[0].set_ylabel("load voltage (V)"); ax[0].legend(ncol=3, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.15))
ax[1].plot(t[m], Y[m, 0] - X[m, 3], color=OI["sky"], lw=0.8, label="measurement error")
ax[1].plot(t[m], XC[m, 3] - X[m, 3], color=OI["orange"], lw=1.2, label="estimate error")
ax[1].set(xlabel="time (s)", ylabel="error (V)"); ax[1].legend(ncol=2, fontsize=8)
fig.savefig(MEDIA / "e3_estimation.png"); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
t0, r0, X0, *_ = simulate(Lq, noise=False)
t1, r1, X1, *_ = simulate(Lq_alt, noise=False)
mm = t0 <= 1.0
ax[0].plot(t0[mm], r0[mm], "--", color=GREY, label="reference r")
ax[0].plot(t0[mm], X0[mm, 3], color=OI["blue"], label="report weights")
ax[0].plot(t1[mm], X1[mm, 3], color=OI["verm"], label="E3.mlx weights")
ax[0].set(xlabel="time (s)", ylabel="load voltage (V)", title="Noise-free tracking, 5 Hz"); ax[0].legend(fontsize=8)
mm = t <= 2
ax[1].plot(t[mm], r[mm], "--", color=GREY, label="reference r")
ax[1].plot(t[mm], X[mm, 3], color=OI["blue"], lw=1, label="true load voltage")
ax[1].set(xlabel="time (s)", title="With process and measurement noise"); ax[1].legend(fontsize=8)
fig.savefig(MEDIA / "e3_tracking.png"); plt.close(fig)

fig, axs = plt.subplots(2, 2, figsize=(9, 5), sharex=True)
for i, (ax, nm) in enumerate(zip(axs.ravel(), ["i1 (A)", "i2 (A)", "iL (A)", "u_load (V)"])):
    ax.plot(t[m], X[m, i], color=OI["blue"], lw=1, label="true"); ax.plot(t[m], XC[m, i], color=OI["orange"], lw=1, label="estimate (updated)")
    ax.set_ylabel(nm)
axs[0, 0].legend(fontsize=8); axs[1, 0].set_xlabel("time (s)"); axs[1, 1].set_xlabel("time (s)")
fig.savefig(MEDIA / "e3_states.png"); plt.close(fig)

# ---------- hero GIF (opaque light card) ----------
mc = {"figure.facecolor": "#f6f7f9", "axes.facecolor": "#f6f7f9", "savefig.facecolor": "#f6f7f9", "savefig.transparent": False,
      "text.color": "#333", "axes.labelcolor": "#333", "xtick.color": "#333", "ytick.color": "#333",
      "axes.edgecolor": "#333", "legend.labelcolor": "#333"}
with mpl.rc_context(mc):
    step, T_GIF = 8, 4.0
    idx = np.arange(0, int(T_GIF / TS), step)
    fig, ax = plt.subplots(2, 1, figsize=(8, 4.6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.2]))
    ax[0].set(xlim=(0, T_GIF), ylim=(-2300, 2300), ylabel="load voltage (V)")
    ax[1].set(ylim=(-900, 900), xlabel="time (s)", ylabel="error (V)")
    ax[0].set_title("LQG on the AC power network: Kalman estimate vs noisy sensor", fontsize=11)
    l_r, = ax[0].plot([], [], "--", color="#555", lw=1.2, label="reference (1 kV, 5 Hz)")
    l_y, = ax[0].plot([], [], color=OI["sky"], lw=0.7, label="noisy measurement")
    l_x, = ax[0].plot([], [], color=OI["blue"], lw=1.6, label="true load voltage")
    l_e, = ax[0].plot([], [], color=OI["verm"], lw=1.6, label="Kalman estimate")
    ax[0].legend(ncol=2, fontsize=8, loc="lower left")
    e_y, = ax[1].plot([], [], color=OI["sky"], lw=0.7, label="measurement error")
    e_x, = ax[1].plot([], [], color=OI["verm"], lw=1.2, label="estimate error")
    ax[1].legend(ncol=2, fontsize=8, loc="lower left")
    for a in ax: a.grid(alpha=0.25)

    def upd(f):
        j = idx[f] + 1
        l_r.set_data(t[:j], r[:j]); l_y.set_data(t[:j], Y[:j, 0]); l_x.set_data(t[:j], X[:j, 3]); l_e.set_data(t[:j], XC[:j, 3])
        e_y.set_data(t[:j], Y[:j, 0] - X[:j, 3]); e_x.set_data(t[:j], XC[:j, 3] - X[:j, 3])
        return l_r, l_y, l_x, l_e, e_y, e_x
    FuncAnimation(fig, upd, frames=len(idx), blit=True).save(MEDIA / "hero_lqg.gif", writer=PillowWriter(fps=15), dpi=70)
    plt.close(fig)
print(cmp_); print(res)
