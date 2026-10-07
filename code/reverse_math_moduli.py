"""
Computable moduli of convergence versus uncomputable limits: the empirical boundary in the series.

Each panel shows a limit used in the series and the error after n finite stages.
Part A  Exact 2D Ising free energy per site on an L x L torus (Kaufman) against the Onsager limit, at K_c and below:
        explicit power-law (critical) or exponential (off-critical) convergence.
Part B  Mass growth of the critical percolation cluster on the diamond lattice (Paper 7): a_n / (p_c Lambda^n) converges
        geometrically, with ratio Lambda_2/Lambda_1 given in closed form.
Part C  van Dam-Hayden embezzlement (Paper 8): 1 - F_n against the explicit asymptotic law (ln 2 - sigma)/H_n.
Part D  Lower approximations Omega_t of the binary-lambda-calculus halting probability (Papers 3, 6).  The certified
        interval [Omega_t, upper bound] does not shrink: no computable modulus is available.

Usage:  python reverse_math_moduli.py [--show]
Figures go to ./figures (PNG, PDF), data to ./data (CSV); override with RM_FIG_DIR, RM_DATA_DIR.
"""
import os
import sys
from functools import lru_cache
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

sys.setrecursionlimit(20000)
FIG_DIR = os.environ.get("RM_FIG_DIR", "figures")
DATA_DIR = os.environ.get("RM_DATA_DIR", "data")
KC = 0.5 * np.log(1 + np.sqrt(2))


def savefig(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=150)


def savedata(name, columns, cols, meta=()):
    os.makedirs(DATA_DIR, exist_ok=True)
    arr = np.column_stack([np.asarray(c, dtype=float) for c in cols])
    with open(os.path.join(DATA_DIR, f"{name}.csv"), "w", encoding="utf-8", newline="\n") as fh:
        for line in meta:
            fh.write(f"# {line}\n")
        fh.write(",".join(columns) + "\n")
        np.savetxt(fh, arr, delimiter=",", fmt="%.12g")


# ================================================================ Part A: Ising free energy
def lnZ_torus(K, L):
    c = np.cosh(2 * K) / np.tanh(2 * K)
    l = np.arange(2 * L)
    gam = np.arccosh(c - np.cos(np.pi * l / L))
    gam[0] = 2 * K + np.log(np.tanh(K))
    logs, signs = [], []
    for idx in (l[1::2], l[0::2]):
        x = L * gam[idx] / 2
        logs.append(np.sum(np.log(2 * np.cosh(x)))); signs.append(1.0)
        sh = 2 * np.sinh(x)
        logs.append(np.sum(np.log(np.abs(sh)))); signs.append(np.prod(np.sign(sh)))
    m = max(logs)
    tot = sum(sg * np.exp(lg - m) for lg, sg in zip(logs, signs))
    return np.log(0.5) + 0.5 * L * L * np.log(2 * np.sinh(2 * K)) + m + np.log(tot)


def onsager(K):
    """Onsager's free energy per site, -beta f = ln(2 cosh 2K) + (1/2pi) int_0^pi ln[(1 + sqrt(1 - k^2 sin^2))/2]."""
    k = 2 * np.sinh(2 * K) / np.cosh(2 * K) ** 2
    val, _ = quad(lambda t: np.log((1 + np.sqrt(max(0.0, 1 - k * k * np.sin(t) ** 2))) / 2), 0, np.pi,
                  limit=400, epsabs=0, epsrel=1e-13)
    return np.log(2 * np.cosh(2 * K)) + val / (2 * np.pi)


def part_a(Ls=(4, 8, 16, 32, 64, 128, 256)):
    out = {}
    for K, lab in ((KC, "K = K_c"), (0.8 * KC, "K = 0.8 K_c")):
        f_inf = onsager(K)
        err = np.array([abs(lnZ_torus(K, L) / L ** 2 - f_inf) for L in Ls])
        out[lab] = err
    Ls = np.array(Ls, float)
    sl = np.polyfit(np.log(Ls[2:]), np.log(out["K = K_c"][2:]), 1)[0]
    rate = -np.polyfit(Ls[:4], np.log(out["K = 0.8 K_c"][:4]), 1)[0]
    # universal amplitude: ln Z_L - L^2 f -> ln Z_CFT(tau = i) = ln[(theta2 + theta3 + theta4) / (2 eta)] for the Ising CFT
    qn = np.exp(-np.pi)
    n = np.arange(0, 40)
    th2 = 2 * np.sum(qn ** ((n + 0.5) ** 2))
    th3 = 1 + 2 * np.sum(qn ** (n[1:] ** 2))
    th4 = 1 + 2 * np.sum((-1.0) ** n[1:] * qn ** (n[1:] ** 2))
    eta = np.exp(-np.pi / 12) * np.prod(1 - np.exp(-2 * np.pi * n[1:]))
    amp_cft = np.log((th2 + th3 + th4) / (2 * eta))
    amp_num = out["K = K_c"][-1] * Ls[-1] ** 2
    print(f"[A] Ising free energy per site, L x L torus vs Onsager: at K_c the error decays as L^{sl:.3f}; "
          f"at 0.8 K_c as exp(-{rate:.3f} L)")
    print(f"    at K_c, L^2 |f_L - f| -> {amp_num:.6f} (L = {int(Ls[-1])}); Ising CFT ln[(th2+th3+th4)/(2 eta)](tau=i) = {amp_cft:.6f}")
    for L, e1, e2 in zip(Ls, out["K = K_c"], out["K = 0.8 K_c"]):
        print(f"    L = {int(L):4d}: |f_L - f| = {e1:.3e} (K_c), {e2:.3e} (0.8 K_c)")
    savedata("ising_free_energy_error", ["L", "err_Kc", "err_0p8Kc"], [Ls, out["K = K_c"], out["K = 0.8 K_c"]],
             ["|ln Z_L / L^2 - Onsager| for the isotropic 2D Ising model on an L x L torus"])
    return Ls, out, sl, rate


# ================================================================ Part B: Paper 7 mass growth
def part_b(levels=40):
    PC = (np.sqrt(5) - 1) / 2
    T = np.array([[12 - 4 * np.sqrt(5), 2 * np.sqrt(5) - 2], [2 * np.sqrt(5) - 4, 2.0]])
    ev = np.sort(np.linalg.eigvals(T).real)[::-1]
    v = np.array([PC, 0.0])
    ratios = []
    for n in range(levels):
        v = T @ v
        ratios.append(v[0] / PC / ev[0] ** (n + 1))
    ratios = np.array(ratios)
    C = ratios[-1]
    err = np.abs(ratios - C)[:-5]
    q = ev[1] / ev[0]
    print(f"[B] diamond-lattice mass recursion: Lambda_1 = {ev[0]:.10f}, Lambda_2 = {ev[1]:.10f}, "
          f"ratio Lambda_2/Lambda_1 = {q:.6f}; |a_n/(p_c Lambda^n) - C| decays by {np.mean(err[1:12] / err[:11]):.6f} per level")
    savedata("dhl_mass_convergence", ["n", "abs_error"], [np.arange(1, len(err) + 1), err],
             [f"T = [[12-4sqrt5, 2sqrt5-2],[2sqrt5-4, 2]]; Lambda2/Lambda1 = {q:.12f}"])
    return err, q


# ================================================================ Part C: embezzlement (Paper 8)
def part_c(kmax=24):
    j = np.arange(1, 2 * 10 ** 7 + 1, dtype=float)
    sigma = np.sum(1 / np.sqrt(2 * j * (2 * j - 1)) - 1 / (2 * j)) + 1 / (8 * j[-1])
    a = np.log(2) - sigma
    ns = 2 ** np.arange(2, kmax + 1)
    errs = []
    for n in ns:
        jj = np.arange(1, n + 1, dtype=float)
        p = 1 / jj; H = p.sum(); p /= H
        q = np.repeat(p / 2, 2)
        F = np.sum(np.sqrt(np.concatenate([p, np.zeros(n)]) * q))
        errs.append(abs((1 - F) - a / H))
    errs = np.array(errs)
    print(f"[C] embezzlement: |(1 - F_n) - (ln2 - sigma)/H_n| from {errs[0]:.2e} (n = 4) to {errs[-1]:.2e} (n = 2^{kmax}); "
          f"times n ln n stays below {np.max(errs * ns * np.log(ns)):.3f}")
    savedata("embezzlement_error", ["n", "abs_error"], [ns, errs], [f"sigma = {sigma:.12f}"])
    return ns, errs


# ================================================================ Part D: Omega approximations
@lru_cache(maxsize=None)
def gen(s, k):
    out = []
    if s >= 2 and s - 2 < k:
        out.append(('v', s - 2))
    if s >= 4:
        out += [('l', b) for b in gen(s - 2, k + 1)]
        for s1 in range(2, s - 3):
            fs = gen(s1, k)
            if not fs:
                continue
            xs = gen(s - 2 - s1, k)
            out += [('a', f, x) for f in fs for x in xs]
    return tuple(out)


def shift(t, d, c=0):
    if t[0] == 'v':
        return ('v', t[1] + d) if t[1] >= c else t
    if t[0] == 'l':
        return ('l', shift(t[1], d, c + 1))
    return ('a', shift(t[1], d, c), shift(t[2], d, c))


def subst(t, j, s):
    if t[0] == 'v':
        return s if t[1] == j else t
    if t[0] == 'l':
        return ('l', subst(t[1], j + 1, shift(s, 1)))
    return ('a', subst(t[1], j, s), subst(t[2], j, s))


def step(t):
    if t[0] == 'a':
        f, x = t[1], t[2]
        if f[0] == 'l':
            return shift(subst(f[1], 0, shift(x, 1)), -1), True
        f2, r = step(f)
        if r:
            return ('a', f2, x), True
        x2, r = step(x)
        return ('a', f, x2), r
    if t[0] == 'l':
        b, r = step(t[1])
        return ('l', b), r
    return t, False


def nodes(t):
    if t[0] == 'v':
        return 1
    if t[0] == 'l':
        return 1 + nodes(t[1])
    return 1 + nodes(t[1]) + nodes(t[2])


def normal_steps(t, budget, max_nodes=3000):
    for k in range(budget + 1):
        t, r = step(t)
        if not r:
            return k
        if nodes(t) > max_nodes:
            return None
    return None


def part_d(tmax=24):
    lengths, steps = [], []
    for L in range(4, tmax + 1):
        for t in gen(L, 0):
            k = normal_steps(t, 16 * tmax)
            lengths.append(L); steps.append(np.inf if k is None else k)
    lengths, steps = np.array(lengths), np.array(steps)
    ts = np.arange(4, tmax + 1)
    lower = np.array([np.sum(2.0 ** -lengths[(lengths <= t) & (steps <= 16 * t)]) for t in ts])
    upper = np.array([lower[i] + np.sum(2.0 ** -lengths[(lengths <= t) & ~(steps <= 16 * t)])
                      + (1.0 - np.sum(2.0 ** -lengths[lengths <= t])) for i, t in enumerate(ts)])
    print(f"[D] Omega_BLC: lower bound rises to {lower[-1]:.6f} at t = {tmax}; certified upper bound {upper[-1]:.6f}; "
          f"certified width {upper[-1] - lower[-1]:.4f} (it was {upper[0] - lower[0]:.4f} at t = 4)")
    savedata("omega_certified_interval", ["t", "lower", "upper"], [ts, lower, upper],
             ["lower: terms of length <= t normalising within 16 t steps; upper: + undecided + Kraft mass of longer terms"])
    return ts, lower, upper


def main():
    plt.rcParams["font.family"] = "DejaVu Sans"
    Ls, iserr, sl, rate = part_a()
    berr, q = part_b()
    ns, eerr = part_c()
    ts, lo, up = part_d()
    fig, ax = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    a = ax[0, 0]
    a.loglog(Ls, iserr["K = K_c"], "o-", label=f"K = K_c (slope {sl:.2f})")
    a.loglog(Ls, iserr["K = 0.8 K_c"], "s-", label="K = 0.8 K_c (exponential)")
    a.set_xlabel("L"); a.set_ylabel("|f_L - f_infinity|"); a.set_title("Ising free energy: explicit modulus")
    a.legend(fontsize=8)
    a = ax[0, 1]
    a.semilogy(np.arange(1, len(berr) + 1), berr, "o-", ms=3, label="|a_n/(p_c Lambda^n) - C|")
    a.semilogy(np.arange(1, len(berr) + 1), berr[0] * q ** np.arange(len(berr)), "--", label=f"(Lambda_2/Lambda_1)^n, ratio {q:.4f}")
    a.set_xlabel("level n"); a.set_title("Diamond-lattice cluster mass (Paper 7): geometric modulus")
    a.legend(fontsize=8)
    a = ax[1, 0]
    a.loglog(ns, eerr, "o-", ms=3, label="|(1 - F_n) - (ln2 - sigma)/H_n|")
    a.loglog(ns, 1 / (ns * np.log(ns)), "--", label="1/(n ln n)")
    a.set_xlabel("n"); a.set_title("Embezzlement (Paper 8): explicit modulus")
    a.legend(fontsize=8)
    a = ax[1, 1]
    a.plot(ts, lo, "o-", label="certified lower bound Omega_t")
    a.plot(ts, up, "s-", label="certified upper bound")
    a.fill_between(ts, lo, up, alpha=0.15)
    a.set_xlabel("stage t"); a.set_ylabel("Omega_BLC"); a.set_title("Halting probability (Papers 3, 6): no computable modulus")
    a.legend(fontsize=8)
    savefig(fig, "reverse_math_moduli")
    print(f"figures written to {FIG_DIR}/, data to {DATA_DIR}/")
    if "--show" in sys.argv:
        plt.show()


if __name__ == "__main__":
    main()
