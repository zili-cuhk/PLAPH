"""Check the TCGA-LUAD preprocessing screen against a working set.

Re-runs the age-, sex- and stage-adjusted death- and censoring-margin Cox screens of
code/R/tcga_preprocessing/run_luad_pilot_dimension_experiment.R (Wald test for each RNA variable;
likelihood-ratio test of a 3-df natural spline for each CNA variable; Efron ties) on the
retained variables of the working set and compares the p-values, coefficients and combined
ranks with those stored in the working-set feature files.

Usage: python3 check_tcga_screen.py WORKING_SET_DIR rna|cna [OUT.csv]
"""
import sys, numpy as np, pandas as pd
from scipy.stats import norm, chi2
W = sys.argv[1]
c = pd.read_csv(W+'/primary_p500_q200_clinical.csv')
Z = pd.read_csv(W+'/primary_p500_q200_Z_RNA_p500.csv'); X = pd.read_csv(W+'/primary_p500_q200_X_CNA_q200.csv')
fz = pd.read_csv(W+'/primary_p500_q200_Z_RNA_p500_features.csv'); fx = pd.read_csv(W+'/primary_p500_q200_X_CNA_q200_features.csv')

def stage_group(s):
    out = np.array([None]*len(s), dtype=object)
    for lab in ['I','II','III','IV']:
        m = s.fillna('').str.contains('STAGE '+lab, regex=False).values
        out[m] = lab
    return out
sg = stage_group(c.stage)
ok = np.isfinite(c.os_months) & (c.os_months > 0) & np.isfinite(c.age) & c.sex.notna() & pd.notna(sg)
ok = ok.values
time = c.os_months.values[ok]; ev = c.os_event.values[ok].astype(float); cen = 1 - ev
cov = np.column_stack([c.age.values[ok], (c.sex.values[ok] == 'Male').astype(float)] +
                      [(sg[ok] == l).astype(float) for l in ['II','III','IV']])
print('screen n', ok.sum(), 'deaths', ev.sum())

def cox_efron(Xm, t, d, iters=50):
    n, p = Xm.shape; b = np.zeros(p)
    order = np.argsort(-t, kind='mergesort'); t, d, Xm = t[order], d[order], Xm[order]
    ut = np.unique(t[d == 1])
    def stats(b):
        eta = Xm @ b; r = np.exp(eta - eta.max())
        ll = 0; g = np.zeros(p); H = np.zeros((p, p))
        for tt in ut:
            risk = t >= tt; D = (t == tt) & (d == 1); m = D.sum()
            s0 = r[risk].sum(); s1 = r[risk] @ Xm[risk]; s2 = (Xm[risk] * r[risk][:, None]).T @ Xm[risk]
            d0 = r[D].sum(); d1 = r[D] @ Xm[D]; d2 = (Xm[D] * r[D][:, None]).T @ Xm[D]
            ll += eta[D].sum()
            g += Xm[D].sum(0)
            for k in range(m):
                f = k / m
                a0 = s0 - f*d0; a1 = s1 - f*d1; a2 = s2 - f*d2
                ll -= np.log(a0) + eta.max()
                g -= a1 / a0
                H -= a2 / a0 - np.outer(a1, a1) / a0**2
        return ll, g, H
    ll, g, H = stats(b)
    for _ in range(iters):
        step = np.linalg.solve(-H, g); bn = b + step
        lln, gn, Hn = stats(bn)
        while lln < ll - 1e-10:
            step /= 2; bn = b + step; lln, gn, Hn = stats(bn)
        conv = abs(lln - ll) < 1e-10
        b, ll, g, H = bn, lln, gn, Hn
        if conv: break
    se = np.sqrt(np.diag(np.linalg.inv(-H)))
    return b, se, ll

def ns3(xv):
    k1, k2 = np.quantile(xv, [1/3, 2/3]); b0, b3 = xv.min(), xv.max()
    knots = [b0, k1, k2, b3]; K = 4
    dk = lambda k: (np.clip(xv - knots[k], 0, None)**3 - np.clip(xv - knots[K-1], 0, None)**3) / (knots[K-1] - knots[k])
    return np.column_stack([xv, dk(0) - dk(2), dk(1) - dk(2)])

which = sys.argv[2] if len(sys.argv) > 2 else 'rna'
rows = []
if which == 'rna':
    for f in fz.feature:
        xv = Z[f].values[ok]; r = {'feature': f}
        for lab, d in [('event', ev), ('censor', cen)]:
            b, se, _ = cox_efron(np.column_stack([xv, cov]), time, d)
            r['p_'+lab] = 2*norm.sf(abs(b[0]/se[0])); r['beta_'+lab] = b[0]
        rows.append(r)
    R = pd.DataFrame(rows).merge(fz, on='feature', suffixes=('_re', ''))
    for col in ['p_event', 'p_censor', 'beta_event', 'beta_censor']:
        rel = np.abs(R[col+'_re'] - R[col]) / np.maximum(np.abs(R[col]), 1e-300)
        print(col, 'max rel diff', rel.max(), 'median', np.median(rel))
else:
    ll0 = {lab: cox_efron(cov, time, d)[2] for lab, d in [('event', ev), ('censor', cen)]}
    for f in fx.feature:
        xv = X[f].values[ok]; r = {'feature': f}
        for lab, d in [('event', ev), ('censor', cen)]:
            _, _, l1 = cox_efron(np.column_stack([xv, cov]), time, d)
            _, _, l3 = cox_efron(np.column_stack([ns3(xv), cov]), time, d)
            r[f'p_{lab}_overall'] = chi2.sf(2*(l3 - ll0[lab]), 3); r[f'p_{lab}_nonlinear'] = chi2.sf(2*(l3 - l1), 2)
        rows.append(r)
    R = pd.DataFrame(rows).merge(fx, on='feature', suffixes=('_re', ''))
    for col in ['p_event_overall', 'p_censor_overall', 'p_event_nonlinear', 'p_censor_nonlinear']:
        rel = np.abs(R[col+'_re'] - R[col]) / np.maximum(np.abs(R[col]), 1e-300)
        print(col, 'max rel diff', rel.max(), 'median', np.median(rel))
R.to_csv(sys.argv[3] if len(sys.argv) > 3 else f'recheck_{which}.csv', index=False)
sc = ('p_event', 'p_censor') if which == 'rna' else ('p_event_overall', 'p_censor_overall')
score = np.maximum(-np.log10(R[sc[0] + '_re']), -np.log10(R[sc[1] + '_re'])).values
same = (np.argsort(-score, kind='stable') == np.argsort(R.combined_rank.values, kind='stable')).all()
print('combined ranking reproduced:', bool(same))
