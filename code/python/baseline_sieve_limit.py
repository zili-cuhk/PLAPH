"""Large-sample limits of the failure coefficients under a piecewise-constant baseline with a fixed
number of intervals.

The simulations generate failure times with baseline hazard h0T(t) = 2t, which a step function with
a fixed number of intervals does not represent exactly. Under independence (tau = 0) the joint
likelihood factorizes, so the failure margin of PLAPH is a Cox model with a piecewise-constant
baseline. For D1 (censoring 0.30) and D4 (censoring 0.50) this script draws one sample of size
200,000, places the interval boundaries at equally spaced quantiles of the observed times (equal
counts, as in the implementation), profiles out the step heights, and reports the ratio of the
fitted failure coefficients of Z1, Z3, Z6 to their true values (1.5, -1.2, 0.9) for 10, 20 and 40
intervals, together with the partial-likelihood fit for reference.

Usage: python3 baseline_sieve_limit.py [OUT.csv] [N]
"""
import csv
import sys
import numpy as np
from scipy.optimize import minimize
from pl_limit_dependent_censoring import simulate, spline_block, cox_fit

TRUE = np.array([1.5, -1.2, 0.9])
IDX = [0, 2, 5]


def step_baseline_fit(Xm, Y, d, k):
    edges = np.quantile(Y, np.linspace(0, 1, k + 1))
    edges[0], edges[-1] = 0.0, np.inf
    expo = np.clip(np.minimum(Y[:, None], edges[None, 1:]) - edges[None, :-1], 0, None)
    bin_of = np.clip(np.searchsorted(edges, Y, side="left") - 1, 0, k - 1)
    events = np.bincount(bin_of, weights=d, minlength=k)

    def negloglik(b):        # baseline heights profiled out: a_k = events_k / sum_i expo_ik exp(eta_i)
        eta = Xm @ b
        r = np.exp(eta)
        s0 = expo.T @ r
        s1 = expo.T @ (r[:, None] * Xm)
        return -(d @ eta - events @ np.log(s0)), -(d @ Xm - (events / s0) @ s1)

    return minimize(negloglik, np.zeros(Xm.shape[1]), jac=True, method="BFGS",
                    options={"gtol": 1e-8, "maxiter": 3000}).x


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "baseline_sieve_limits.csv"
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    rows = []
    for label, cs, seed in [("D1 (cens 0.30, tau 0)", 15.9, 4100), ("D4 (cens 0.50, tau 0)", 1.3, 4103)]:
        Y, d, Z, X = simulate(N, 0.0, cs, 1.5, seed=seed)
        Xm = np.column_stack([Z] + [spline_block(X[:, j]) for j in range(4)])
        row = {"design": label, "censoring_rate": 1 - d.mean()}
        for k in (10, 20, 40):
            b = step_baseline_fit(Xm, Y, d, k)
            row[f"ratio_k{k}"] = float(np.mean(b[IDX] / TRUE))
        b, _ = cox_fit(Xm, Y, d)
        row["ratio_partial_likelihood"] = float(np.mean(b[IDX] / TRUE))
        rows.append(row)
        print(label, {k: round(v, 3) for k, v in row.items() if k != "design"}, flush=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print("written", out)
