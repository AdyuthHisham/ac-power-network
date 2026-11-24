"""Shared model, design routines and plot style.

Python re-implementation of the MATLAB model from the SSY285 assignments (E1-E3).
Parameters are the ones used in the MATLAB scripts (R1=30, L1=0.1, C1=0.5 uF, ...).
"""
import json
import pathlib

import matplotlib as mpl
import numpy as np
from scipy.linalg import expm, solve_discrete_are

ROOT = pathlib.Path(__file__).resolve().parents[1]
MEDIA = ROOT / "media"
METRICS = ROOT / "docs" / "metrics"
MEDIA.mkdir(exist_ok=True)
METRICS.mkdir(parents=True, exist_ok=True)

# ---- plant (E1) ----
R1, L1, C1 = 30.0, 0.1, 0.5e-6      # line 1: 0.3 ohm/km, 1 mH/km, 100 km, C split over two ends
R2, L2, C2 = 30.0, 0.2, 1.5e-6      # line 2: 0.15 ohm/km, 1 mH/km, 200 km
RL, LL = 100.0, 1.0                 # series RL load
E, F = 10_000.0, 50.0
W = 2 * np.pi * F
TS = 0.005


def plant(rl=RL, ll=LL):
    """x = [i1, i2, iL, uL], u = [u1, u2], y = uL."""
    cs = C1 + C2
    A = np.array([[-R1 / L1, 0, 0, -1 / L1],
                  [0, -R2 / L2, 0, -1 / L2],
                  [0, 0, -rl / ll, 1 / ll],
                  [1 / cs, 1 / cs, -1 / cs, 0]])
    B = np.array([[1 / L1, 0], [0, 1 / L2], [0, 0], [0, 0]])
    C = np.array([[0.0, 0, 0, 1]])
    return A, B, C


def zoh(A, B, ts=TS):
    n, m = B.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n], M[:n, n:] = A, B
    Md = expm(M * ts)
    return Md[:n, :n], Md[:n, n:]


def dlqr(Ad, Bd, Qx, Qu):
    P = solve_discrete_are(Ad, Bd, Qx, Qu)
    return np.linalg.solve(Qu + Bd.T @ P @ Bd, Bd.T @ P @ Ad)


def kalman_gain(Ad, Bd, Cd, R1, R2):
    """Predictor-form steady-state gain (x+ = Ad x + Bd u + K (y - Cd x))."""
    Q = Bd @ R1 @ Bd.T            # Nd = Bd
    P = solve_discrete_are(Ad.T, Cd.T, Q, R2)
    K = Ad @ P @ Cd.T @ np.linalg.inv(Cd @ P @ Cd.T + R2)
    return K, P


def _enc(o):
    if isinstance(o, (complex, np.complexfloating)):
        return {"re": float(o.real), "im": float(o.imag)}
    if isinstance(o, np.ndarray):
        return o.tolist()
    return float(o)


def save_metrics(name, d):
    (METRICS / f"{name}.json").write_text(json.dumps(d, indent=2, default=_enc))


# ---- plot style: colour-blind-safe (Okabe-Ito), mid-grey text, transparent background ----
OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
          sky="#56B4E9", pink="#CC79A7")
GREY = "#808080"


def style():
    mpl.rcParams.update({
        "figure.facecolor": "none", "axes.facecolor": "none", "savefig.facecolor": "none",
        "savefig.transparent": True, "text.color": GREY, "axes.labelcolor": GREY,
        "axes.edgecolor": GREY, "xtick.color": GREY, "ytick.color": GREY,
        "grid.color": GREY, "grid.alpha": 0.25, "axes.grid": True, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False, "legend.labelcolor": GREY,
        "axes.prop_cycle": mpl.cycler(color=[OI["blue"], OI["orange"], OI["green"], OI["verm"], OI["sky"]]),
        "font.size": 10, "savefig.dpi": 140, "savefig.bbox": "tight",
    })
