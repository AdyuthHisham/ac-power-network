"""E1: sinusoidal steady state, step response, initial-condition response, poles and zeros."""
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from common import *

style()
A, B, C = plant()
sys = signal.StateSpace(A, B, C, np.zeros((1, 2)))
U = lambda t: np.column_stack([np.sqrt(2) * E * np.sin(W * t), np.sqrt(2) * E * np.cos(W * t)])

# --- steady state: analytic phasor vs time simulation ---
G = lambda s: (C @ np.linalg.solve(s * np.eye(4) - A, B))[0]
g1, g2 = G(1j * W)
phasor = np.sqrt(2) * E * (-1j * g1 + g2)          # y_ss(t) = Re(phasor e^{jwt})
t = np.linspace(0, 0.3, 6001)
_, y, _ = signal.lsim(sys, U(t), t)
tail = t > 0.2
k = np.argmax(np.abs(y[tail]))
sim_amp = np.max(np.abs(y[tail]))

fig, ax = plt.subplots(figsize=(7, 3.2))
ax.plot(t * 1e3, y / 1e3, label="load voltage (simulation)")
ax.plot(t * 1e3, np.abs(phasor) / 1e3 * np.cos(W * t + np.angle(phasor)), "--", color=OI["orange"], label="analytic steady state")
ax.set(xlabel="time (ms)", ylabel="u_load (kV)", title="Sinusoidal steady state of the load voltage")
ax.legend(); fig.savefig(MEDIA / "e1_steady_state.png"); plt.close(fig)

# --- linearised model: Jacobians of a linear system are A, B, C (checked by finite differences) ---
f = lambda x, u: A @ x + B @ u
eps = 1e-3
JA = np.column_stack([(f(eps * e, np.zeros(2)) - f(-eps * e, np.zeros(2))) / (2 * eps) for e in np.eye(4)])
jac_err = float(np.max(np.abs(JA - A)))

t = np.linspace(0, 0.1, 2001)
_, ys, _ = signal.lsim(sys, np.ones((t.size, 2)), t)
fig, ax = plt.subplots(figsize=(7, 3.2))
ax.plot(t * 1e3, ys / 1e3)
ax.set(xlabel="time (ms)", ylabel="u_load (kV)", title="Linearised model, unit step on u1 and u2")
fig.savefig(MEDIA / "e1_step.png"); plt.close(fig)

# --- initial load voltage, 5 Hz / 1 kV input ---
t = np.linspace(0, 0.25, 5001)
d = 1000 * np.column_stack([np.sin(2 * np.pi * 5 * t), np.cos(2 * np.pi * 5 * t)])
fig, ax = plt.subplots(figsize=(7, 3.4))
for v0 in (-1000, -750, 0, 250, 1000):
    _, yy, _ = signal.lsim(sys, d, t, X0=[0, 0, 0, v0])
    ax.plot(t * 1e3, yy, label=f"u_load(0) = {v0} V")
ax.set(xlabel="time (ms)", ylabel="u_load (V)", title="Response to a 5 Hz, 1 kV input for five initial load voltages")
ax.legend(fontsize=8, ncol=2); fig.savefig(MEDIA / "e1_initial.png"); plt.close(fig)

# --- transfer functions, poles, zeros ---
poles = np.linalg.eigvals(A)
tf = [signal.ss2tf(A, B, C, np.zeros((1, 2)), input=i) for i in (0, 1)]
clean = lambda n: np.where(np.abs(n) < 1e-9 * np.abs(n).max(), 0.0, n)   # drop numerical-zero leading terms
tf = [(clean(num[0])[None, :], den) for num, den in tf]
zeros = [np.roots(np.trim_zeros(num[0], "f")) for num, _ in tf]
fig, ax = plt.subplots(figsize=(5.2, 3.6))
ax.plot(poles.real, poles.imag, "x", ms=9, mew=2, label="poles")
for i, z in enumerate(zeros):
    ax.plot(z.real, z.imag, "o", mfc="none", ms=9 + 4 * i, mew=1.8, label=f"zeros, input {i + 1}")
ax.set(xlabel="real (1/s)", ylabel="imag (1/s)", title="Poles and transmission zeros")
ax.legend(fontsize=8); fig.savefig(MEDIA / "e1_pzmap.png"); plt.close(fig)

save_metrics("e1", dict(
    steady_state_amplitude_V=abs(phasor), steady_state_phase_deg=np.degrees(np.angle(phasor)),
    steady_state_sim_peak_V=float(sim_amp),
    poles=[complex(p) for p in poles], zeros_input1=[complex(z) for z in zeros[0]],
    zeros_input2=[complex(z) for z in zeros[1]],
    tf_input1_num=tf[0][0][0].tolist(), tf_input2_num=tf[1][0][0].tolist(), tf_den=tf[0][1].tolist(),
    jacobian_max_abs_error=jac_err, max_pole_real=float(poles.real.max())))
print("phasor", phasor, abs(phasor), np.degrees(np.angle(phasor)), "sim", sim_amp)
print("poles", poles); print("zeros", zeros); print(tf)
