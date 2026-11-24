"""E2: controllability / observability for symbolic RL, LL; detectability; conditioning; ZOH discretisation."""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

from common import *

style()
# ---------- symbolic part ----------
rl, ll = sp.symbols("RL LL", positive=True)
cs = sp.Rational(1, 2) / 10**6 + sp.Rational(3, 2) / 10**6
As = sp.Matrix([[-300, 0, 0, -10], [0, -150, 0, -5], [0, 0, -rl / ll, 1 / ll], [1 / cs, 1 / cs, -1 / cs, 0]])
Bs = sp.Matrix([[10, 0], [0, 5], [0, 0], [0, 0]])
Cs = sp.Matrix([[0, 0, 0, 1]])
Obs = sp.Matrix.vstack(*[Cs * As**k for k in range(4)])
det_obs = sp.factor(sp.simplify(Obs.det()))
Ctr = sp.Matrix.hstack(*[As**k * Bs for k in range(4)])
gram_det = sp.factor(sp.simplify((Ctr * Ctr.T).det()))
roots = sp.solve(sp.Eq(det_obs, 0), rl)
print("det(Obs) =", det_obs); print("unobservable at RL =", roots); print("det(Ctr Ctr^T) =", gram_det)

# ---------- numeric PBH scan over RL/LL ----------
ratios = np.logspace(1, 3, 4001)


def pbh_min_sigma(A, C, extra=None):
    lam = np.linalg.eigvals(A)
    return min(np.linalg.svd(np.vstack([A - l * np.eye(4), C]), compute_uv=False)[-1] for l in lam)


obs_sig = []
ctr_sig = []
for r in ratios:
    A, B, C = plant(rl=r, ll=1.0)
    obs_sig.append(pbh_min_sigma(A, C))
    lam = np.linalg.eigvals(A)
    ctr_sig.append(min(np.linalg.svd(np.hstack([A - l * np.eye(4), B]), compute_uv=False)[-1] for l in lam))
obs_sig, ctr_sig = np.array(obs_sig), np.array(ctr_sig)

# rank checks at the two bad ratios and at the nominal load
rank_info = {}
for r in (100.0, 150.0, 300.0):
    A, B, C = plant(rl=r, ll=1.0)
    O = np.vstack([C @ np.linalg.matrix_power(A, k) for k in range(4)])
    rank_info[f"RL={r:g},LL=1"] = dict(
        pbh_obs_min_sigma=pbh_min_sigma(A, C),
        eig_A=np.sort(np.linalg.eigvals(A).real).tolist())

# detectability: K = [0 0 1/LL 0]' removes the coupling of x3 so the unobservable mode is the stable load pole
det_ok = {}
for r in (150.0, 300.0):
    A, B, C = plant(rl=r, ll=1.0)
    K = np.array([[0], [0], [1.0 / 1.0], [0]])
    det_ok[f"RL={r:g}"] = dict(eig_A_minus_KC=np.sort(np.linalg.eigvals(A - K @ C).real).tolist(),
                               stable=bool(np.all(np.linalg.eigvals(A - K @ C).real < 0)))

fig, ax = plt.subplots(figsize=(7, 3.4))
ax.loglog(ratios, obs_sig, label="observability (PBH, min singular value)")
ax.loglog(ratios, ctr_sig, label="controllability (PBH, min singular value)")
for r_, c_ in ((150, OI["orange"]), (300, OI["verm"])):
    ax.axvline(r_, color=c_, ls="--", lw=1)
    ax.text(r_ * 1.03, ax.get_ylim()[0] * 5, f"RL/LL = {r_}", color=c_, rotation=90, va="bottom", fontsize=8)
ax.set(xlabel="RL / LL (1/s)", ylabel="min over modes of sigma_min", title="Where the load makes the output unobservable")
ax.legend(fontsize=8); fig.savefig(MEDIA / "e2_pbh.png"); plt.close(fig)

# ---------- conditioning and discretisation ----------
A, B, C = plant(); Ad, Bd = zoh(A, B)
ctr = lambda a, b: np.hstack([np.linalg.matrix_power(a, k) @ b for k in range(4)])
obs = lambda a, c: np.vstack([c @ np.linalg.matrix_power(a, k) for k in range(4)])
cond = dict(ctrb_continuous=np.linalg.cond(ctr(A, B)), obsv_continuous=np.linalg.cond(obs(A, C)),
            ctrb_discrete=np.linalg.cond(ctr(Ad, Bd)), obsv_discrete=np.linalg.cond(obs(Ad, C)))
lam_c, lam_d = np.linalg.eigvals(A), np.linalg.eigvals(Ad)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw=dict(width_ratios=[1, 1.1]))
th = np.linspace(0, 2 * np.pi, 400)
a1.plot(np.cos(th), np.sin(th), color=GREY, lw=1)
a1.plot(lam_d.real, lam_d.imag, "x", ms=10, mew=2, color=OI["blue"])
a1.set(aspect="equal", xlabel="real", ylabel="imag", title=f"Discrete eigenvalues (Ts = {TS * 1e3:g} ms)")
names = ["ctrb\ncont.", "obsv\ncont.", "ctrb\ndisc.", "obsv\ndisc."]
vals = list(cond.values())
a2.bar(names, vals, color=[OI["blue"], OI["orange"], OI["blue"], OI["orange"]])
a2.set_yscale("log"); a2.set(ylabel="condition number", title="Conditioning of the test matrices")
for i, v in enumerate(vals):
    a2.text(i, v * 1.3, f"{v:.3g}", ha="center", fontsize=8)
a2.set_ylim(top=max(vals) * 40)
fig.savefig(MEDIA / "e2_discrete.png"); plt.close(fig)

save_metrics("e2", dict(
    det_observability_matrix=str(det_obs), unobservable_RL_in_terms_of_LL=[str(r) for r in roots],
    det_ctrb_gram=str(gram_det), rank_checks=rank_info, detectability_check=det_ok,
    cond=cond, Ad=Ad, Bd=Bd, discrete_eig=lam_d, discrete_eig_abs=np.abs(lam_d)))
print(cond); print(rank_info); print(det_ok)
