"""Large-sample limits of the failure-margin partial-likelihood (PL) fit under dependent censoring.

For the simulation designs of the article, draws one large sample (default n = 200,000) from the
data-generating process of plaph_run_one_scenario_impl.m (discriminating design,
"shared_balanced_bounded" additive functions, Frank copula) and fits the failure-margin Cox model
that treats censoring as independent: linear terms for Z1..Z10 and centered quadratic B-splines
(four interior quantile knots) for X1..X4. The fitted coefficients approximate the PL limits.
Under tau > 0 the censoring-only covariates Z2, Z4, Z7 have nonzero limits; the script reports
them, the limits of the failure coefficients (true values 1.5, -1.2, 0.9), and the implied ratio
limit / standard error at the sample sizes of the main designs (n = 500, 700).

Usage: python3 pl_limit_dependent_censoring.py [OUT.csv] [N]
"""
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm
from scipy.interpolate import BSpline


def frank_theta(tau):
    if abs(tau) < 1e-12:
        return 0.0
    if tau < 0:
        return -frank_theta(-tau)
    d1 = lambda a: quad(lambda t: t / np.expm1(t), 0, a)[0] / a
    return brentq(lambda a: 1 - 4 / a * (1 - d1(a)) - tau, 1e-6, 300)


def frank_draw(theta, n, rng):
    u, w = rng.random(n), rng.random(n)
    if theta == 0:
        return u, w
    v = -np.log1p(w * np.expm1(-theta) / (w + (1 - w) * np.exp(-theta * u))) / theta
    return u, v


def additive(X):
    """Failure and censoring additive predictors of the shared_balanced_bounded design (scale 2)."""
    f = np.column_stack([0.8 * X[:, 0], np.sqrt(5 / 2) * (6 * X[:, 1] ** 2 - 6 * X[:, 1] + 1),
                         np.sin(2 * np.pi * X[:, 2]), np.cos(2 * np.pi * X[:, 3])]) * 2.0
    return f.sum(1), (f * np.array([-1, 1, -1, -1])).sum(1)


def simulate(n, tau, cs, gscale, seed, p=10, rho=0.25):
    rng = np.random.default_rng(seed)
    idx = np.arange(p)
    Z = rng.multivariate_normal(np.zeros(p), rho ** np.abs(idx[:, None] - idx[None, :]), n)
    q = 6
    jdx = np.arange(q)
    L = rng.multivariate_normal(np.zeros(q), rho ** np.abs(jdx[:, None] - jdx[None, :]), n)
    X = L.copy()
    X[:, 1:6] = norm.cdf(L[:, 1:6])
    beta = np.zeros(p); gamma = np.zeros(p)
    beta[[0, 2, 5]] = 1.5 * np.array([1, -0.8, 0.6])
    gamma[[1, 3, 6]] = gscale * np.array([-0.8, 1, 0.6])
    aT, aC = additive(X)
    u, v = frank_draw(frank_theta(tau), n, rng)
    T = np.sqrt(-np.log(1 - u) / np.exp(Z @ beta + aT))       # h0T(t) = 2t
    C = -cs * np.log(1 - v) / np.exp(Z @ gamma + aC)           # exponential censoring
    return np.minimum(T, C), (T <= C).astype(float), Z, X


def spline_block(x, interior=4, order=3):
    k = order - 1
    t = np.quantile(x, np.linspace(0, 1, interior + 2))
    knots = np.r_[[t[0]] * order, t[1:-1], [t[-1]] * order]
    B = BSpline.design_matrix(np.clip(x, t[0], t[-1]), knots, k).toarray()
    B = B - B.mean(0)
    return B[:, :-1]                                           # centered basis has rank K-1


def cox_fit(Xm, Y, d, iters=50):
    """Breslow partial-likelihood Newton fit with O(n p^2) memory."""
    o = np.argsort(-Y, kind="mergesort")
    Xm, d = Xm[o], d[o]
    b = np.zeros(Xm.shape[1])
    for _ in range(iters):
        eta = Xm @ b
        r = np.exp(eta - eta.max())
        s0 = np.cumsum(r)
        s1 = np.cumsum(r[:, None] * Xm, 0)
        ev = d == 1
        mu = s1[ev] / s0[ev, None]
        g = (Xm[ev] - mu).sum(0)
        # sum over events of S2/S0 equals X' diag(r * c) X, c_j = sum_{events i, i >= j in risk order} 1/S0_i
        inc = np.where(ev, 1.0 / s0, 0.0)
        c = np.cumsum(inc[::-1])[::-1]
        H = (Xm * (r * c)[:, None]).T @ Xm - mu.T @ mu
        step = np.linalg.solve(H, g)
        b += step
        if np.abs(step).max() < 1e-9:
            break
    return b, np.sqrt(np.diag(np.linalg.inv(H)))


DESIGNS = [  # label, tau, censoring scale, censoring-only coefficient scale
    ("D1 (cens 0.30, tau 0)", 0.0, 15.9, 1.5), ("D2 (cens 0.30, tau 0.5)", 0.5, 15.7, 1.5),
    ("D3 (cens 0.30, tau 0.8)", 0.8, 15.7, 1.5), ("D4 (cens 0.50, tau 0)", 0.0, 1.3, 1.5),
    ("D5 (cens 0.50, tau 0.5)", 0.5, 1.3, 1.5), ("D6 (cens 0.50, tau 0.8)", 0.8, 1.3, 1.5)]
NS = (500, 700)                                                # sample sizes of the main designs

if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "pl_limits.csv"
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    rows = []
    for k, (label, tau, cs, gs) in enumerate(DESIGNS):
        Y, d, Z, X = simulate(N, tau, cs, gs, seed=4100 + k)
        Xm = np.column_stack([Z] + [spline_block(X[:, j]) for j in range(4)])
        b, se = cox_fit(Xm, Y, d)
        unit = se * np.sqrt(N)                                 # standard error at n = 1
        row = {"design": label, "tau": tau, "censoring_rate": 1 - d.mean(), "censoring_only_scale": gs}
        for j, name in [(1, "Z2"), (3, "Z4"), (6, "Z7")]:
            row[f"limit_{name}"] = b[j]
            for n in NS:
                row[f"z_{name}_n{n}"] = abs(b[j]) / unit[j] * np.sqrt(n)
        for j, name, true in [(0, "Z1", 1.5), (2, "Z3", -1.2), (5, "Z6", 0.9)]:
            row[f"limit_{name}"] = b[j]
            row[f"ratio_{name}"] = b[j] / true
        rows.append(row)
        print(f"{label:30s} cens={row['censoring_rate']:.3f}  Z2,Z4,Z7 = "
              f"{b[1]:+.3f} {b[3]:+.3f} {b[6]:+.3f}   failure ratios = "
              f"{row['ratio_Z1']:.3f} {row['ratio_Z3']:.3f} {row['ratio_Z6']:.3f}   "
              f"max z at n=500/700 = " +
              "/".join(f"{max(row[f'z_{m}_n{n}'] for m in ('Z2','Z4','Z7')):.1f}" for n in NS))
    import csv
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print("written", out)
