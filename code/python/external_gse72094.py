"""Prespecified external evaluation of the TCGA failure-margin RNA effect in GSE72094.

Score R_i = sum_g beta_g^TCGA * Xbar_ig, where each probe of a selected gene is
standardized across external samples and probes of one gene are averaged.
Cox model (Efron ties) for OS on the score, unadjusted and adjusted for age,
male sex and stage II/III/IV (stage I reference); Harrell's C-index.
"""
import sys, json
import numpy as np, pandas as pd
from math import erf, sqrt

BETA_TCGA = {'CYP17A1': -0.395}            # tau = 0 primary fit (gamma = 1)
PROBES = {'CYP17A1': ['probe_merck-NM_000102_at']}


def cox_efron(time, event, X, tol=1e-10, maxit=50):
    X = np.asarray(X, float); n, p = X.shape
    order = np.argsort(-time, kind='mergesort'); t, d, Z = time[order], event[order], X[order]
    b = np.zeros(p)
    for _ in range(maxit):
        eta = Z @ b; w = np.exp(eta - eta.max())
        S0 = np.cumsum(w); S1 = np.cumsum(w[:, None] * Z, 0); S2 = np.cumsum(w[:, None, None] * Z[:, :, None] * Z[:, None, :], 0)
        ll, U, H = 0.0, np.zeros(p), np.zeros((p, p))
        i = 0
        while i < n:
            j = i
            while j + 1 < n and t[j + 1] == t[i]: j += 1
            last = j                                 # risk set = rows 0..last (descending time)
            ev = [k for k in range(i, j + 1) if d[k]]
            m = len(ev)
            if m:
                wd = w[ev]; Zd = Z[ev]
                D0, D1 = wd.sum(), (wd[:, None] * Zd).sum(0)
                D2 = (wd[:, None, None] * Zd[:, :, None] * Zd[:, None, :]).sum(0)
                for r in range(m):
                    f = r / m
                    a0 = S0[last] - f * D0; a1 = S1[last] - f * D1; a2 = S2[last] - f * D2
                    ll -= np.log(a0) + eta.max()
                    U -= a1 / a0
                    H -= a2 / a0 - np.outer(a1, a1) / a0 ** 2
                ll += (eta[ev]).sum(); U += Zd.sum(0)
            i = j + 1
        step = np.linalg.solve(-H, U)
        b = b + step
        if np.max(np.abs(step)) < tol: break
    cov = np.linalg.inv(-H)
    return b, np.sqrt(np.diag(cov)), ll


def pnorm2(z):
    return 2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2))))


def harrell_c(time, event, risk):
    num = den = 0.0
    for i in np.where(event == 1)[0]:
        comp = time > time[i]
        den += comp.sum()
        num += (risk[i] > risk[comp]).sum() + 0.5 * (risk[i] == risk[comp]).sum()
    return num / den


def selftest(seed=1):
    rng = np.random.default_rng(seed); n = 4000
    x = rng.normal(size=(n, 2)); t = rng.exponential(1 / np.exp(0.7 * x[:, 0] - 0.4 * x[:, 1]))
    c = rng.exponential(1.5, n); time = np.round(np.minimum(t, c), 2); ev = (t <= c).astype(int)
    b, se, _ = cox_efron(time, ev, x)
    assert abs(b[0] - 0.7) < 4 * se[0] and abs(b[1] + 0.4) < 4 * se[1], (b, se)
    return b, se


if __name__ == '__main__':
    print('self-test (truth 0.7, -0.4):', selftest())
    D = pd.read_csv(sys.argv[1], index_col=0)
    D['time'] = pd.to_numeric(D.survival_time_in_days, errors='coerce') / 30.4375
    D['event'] = D.vital_status.map({'Dead': 1, 'Alive': 0})
    Xg = {}
    for g, pr in PROBES.items():
        Z = (D[pr] - D[pr].mean()) / D[pr].std(ddof=1)
        Xg[g] = Z.mean(axis=1)
    D['score'] = sum(BETA_TCGA[g] * Xg[g] for g in BETA_TCGA)
    E = D[D.time.notna() & D.event.notna() & D.score.notna() & (D.time > 0)].copy()
    t, ev = E.time.to_numpy(), E.event.to_numpy().astype(int)
    b, se, _ = cox_efron(t, ev, E[['score']].to_numpy())
    out = {'n': len(E), 'deaths': int(ev.sum()), 'hr': float(np.exp(b[0])),
           'ci': [float(np.exp(b[0] - 1.96 * se[0])), float(np.exp(b[0] + 1.96 * se[0]))],
           'p': pnorm2(b[0] / se[0]), 'c_index': harrell_c(t, ev, E.score.to_numpy())}
    # gene-level external coefficient (standardized expression), unadjusted
    genes = {}
    for g in BETA_TCGA:
        bg, sg, _ = cox_efron(t, ev, Xg[g].loc[E.index].to_numpy()[:, None])
        genes[g] = {'tcga': BETA_TCGA[g], 'external': float(bg[0]), 'se': float(sg[0]), 'p': pnorm2(bg[0] / sg[0]),
                    'sign_agree': bool(np.sign(bg[0]) == np.sign(BETA_TCGA[g]))}
    # adjusted model
    st = E.Stage.astype(str).str.upper()
    E['age'] = pd.to_numeric(E.age_at_diagnosis, errors='coerce')
    E['male'] = E.gender.map({'M': 1, 'F': 0})
    E['s2'] = st.str.startswith('2').astype(float); E['s3'] = st.str.startswith('3').astype(float); E['s4'] = st.str.startswith('4').astype(float)
    E.loc[~st.str.match(r'^[1-4]'), ['s2', 's3', 's4']] = np.nan
    A = E.dropna(subset=['age', 'male', 's2', 's3', 's4'])
    ba, sa, _ = cox_efron(A.time.to_numpy(), A.event.to_numpy().astype(int), A[['score', 'age', 'male', 's2', 's3', 's4']].to_numpy())
    out.update({'adj_n': len(A), 'adj_deaths': int(A.event.sum()), 'adj_hr': float(np.exp(ba[0])),
                'adj_ci': [float(np.exp(ba[0] - 1.96 * sa[0])), float(np.exp(ba[0] + 1.96 * sa[0]))],
                'adj_p': pnorm2(ba[0] / sa[0]), 'genes': genes,
                'expr_summary': {'median': float(D[PROBES['CYP17A1'][0]].median()), 'iqr': [float(D[PROBES['CYP17A1'][0]].quantile(q)) for q in (.25, .75)]}})
    print(json.dumps(out, indent=1))
    json.dump(out, open(sys.argv[2], 'w'), indent=1)
