"""Rows of Supplementary Table tab:tau-misspec (misspecified working association).

Usage: python3 fill_tau_misspec.py RESULTS_SIMULATION_DIR [SUPP_TEX]
D5 (tau = 0.5) and D6 (tau = 0.8) data, fitted with the Frank copula at the working tau shown;
working tau 0 is PLAPH-ind (ablation_tau0), the true value is PLAPH (main), the others are
tau_misspec. With SUPP_TEX, compares the rows with the file and reports differences.
"""
import glob, os, re, sys
import numpy as np, pandas as pd

res = sys.argv[1]
B = np.array([1.5, -1.2, 0.9]); G = np.array([-1.2, 1.5, 0.9])


def load(pat):
    fs = glob.glob(os.path.join(res, pat, 'chunk_*', 'scenario_001_replicates.csv'))
    return pd.concat([pd.read_csv(f) for f in fs])


def row(label, wt, R):
    ex = lambda p: ((R[p + '_fp'] == 0) & (R[p + '_fn'] == 0))
    all4 = ex('beta') & ex('phi') & ex('mpl_add_T') & ex('mpl_add_C')
    mse = lambda pre, T: np.mean([np.mean((R[f'{pre}{k+1}'].fillna(0) - T[k]) ** 2) for k in range(3)])
    cp = lambda pre: np.mean([np.where(np.isfinite(R[f'{pre}{k+1}']), R[f'{pre}{k+1}'], 0)[all4.values].mean()
                              for k in range(3)])
    v = [f"{R.mpl_accepted.mean():.3f}", f"{all4.mean():.3f}", f"{ex('beta').mean():.3f}", f"{ex('phi').mean():.3f}",
         f"{mse('mpl_beta_active', B):.4f}", f"{mse('mpl_phi_active', G):.4f}",
         f"{R.imse_mpl_add_T.mean():.4f}", f"{R.imse_mpl_add_C.mean():.4f}",
         f"{cp('mpl_beta_cp_active'):.3f}", f"{cp('mpl_phi_cp_active'):.3f}"]
    return f"{label} & {wt} & " + " & ".join(v) + r" \\"


rows = [row('D5', '0', load('ablation_tau0/*_S05_*')), row('', '0.2', load('tau_misspec/*S101*')),
        row('', '0.5', load('main/V2_main_S05_*')), row('', '0.8', load('tau_misspec/*S102*')),
        r'\addlinespace',
        row('D6', '0', load('ablation_tau0/*_S06_*')), row('', '0.5', load('tau_misspec/*S103*')),
        row('', '0.8', load('main/V2_main_S06_*'))]
print('\n'.join(rows))
if len(sys.argv) > 2:
    t = open(sys.argv[2]).read(); a = t.index(r'\label{tab:tau-misspec}')
    m = t.index(r'\midrule', a) + len(r'\midrule'); b = t.index(r'\bottomrule', a)
    have = [L.strip() for L in t[m:b].strip().split('\n') if L.strip()]
    bad = sum(re.sub(r'\s+', '', x) != re.sub(r'\s+', '', y) for x, y in zip(rows, have)) + (len(rows) != len(have))
    print(f'{bad} differences')
