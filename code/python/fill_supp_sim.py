"""Fill simulation tables of PLAPH_supplement.tex from replicate-level CSVs.

Definitions (as stated in the supplement):
* Exact recovery, FP, FN, size: means over all 200 replications.
* Linear coefficients, unconditional: all replications; an omitted coefficient
  is 0 in bias/ESE/MSE and its interval is noncovering; ASE averages finite
  model-based SEs of selected coefficients.
* Conditional (superscript c): replications in which all four supports are
  recovered exactly (PLAPH, PLAPH-ind); for PL and grpsurv, the two
  failure-margin supports.
"""
import glob, re, sys
import numpy as np, pandas as pd

runs, tex = sys.argv[1], sys.argv[2]
BETA = np.array([1.5, -1.2, 0.9]); GAMMA = np.array([-1.2, 1.5, 0.9])


def load(exp, s):
    fs = sorted(glob.glob(f'{runs}/{exp}/V2_{exp}_S{s:02d}_*/chunk_*/scenario_001_replicates.csv'))
    if not fs:
        return None
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True).sort_values('seed')


def load_grp(s):
    fs = sorted(f for f in glob.glob(f'{runs}/grpsurv_results/V2_grpsurv_export_S{s:02d}_*/rows_*.csv')
                if not f.endswith('_curves.csv'))
    if not fs:
        return None
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True).sort_values('seed')


MAIN = {s: load('main', s) for s in range(1, 14)}
IND = {s: load('ablation_tau0', s) for s in range(1, 14)}
GRP = {s: load_grp(s) for s in range(1, 14)}
TAU = {s: t for s, t in zip(range(1, 14), [0, .5, .8, 0, .5, .8, 0, .5, .8, .5, .5, 0, .8])}


def exact(R, a, b):
    return (R[a] == 0) & (R[b] == 0)


def all4(R):
    return exact(R, 'beta_fp', 'beta_fn') & exact(R, 'phi_fp', 'phi_fn') & \
        exact(R, 'mpl_add_T_fp', 'mpl_add_T_fn') & exact(R, 'mpl_add_C_fp', 'mpl_add_C_fn')


f3 = lambda x: f'{x:.3f}'
f4 = lambda x: f'{x:.4f}'


def sci(x):
    if not np.isfinite(x):
        return '--'
    m, e = f'{x:.1e}'.split('e')
    return rf'${m}\times10^{{{int(e)}}}$'


# ---------------------------------------------------------------- diagnostics
def diag_row(s):
    R = MAIN[s]; acc = R.mpl_accepted == 1; A = R[acc]
    agree = [f3((A[c] == 1).mean()) for c in ['mpl_pilot_final_agree_linear_T', 'mpl_pilot_final_agree_linear_C',
                                             'mpl_pilot_final_agree_add_T', 'mpl_pilot_final_agree_add_C']]
    ind = '--' if TAU[s] == 0 else f3((IND[s].mpl_accepted == 1).mean())
    tie = (R.tuning_tie_count == R.tuning_valid_count).mean()
    v = [f3(R.screen_all_sure.mean()), f3(acc.mean()), ind, f3((R.pl_accepted == 1).mean())] + agree + \
        [f3(tie), f3(R.tuning_at_boundary.mean()), sci(A.mpl_refit_stationarity_inf.max())]
    return f'D{s} & ' + ' & '.join(v) + r' \\'


# ---------------------------------------------------------- selection detail
def sel_vals(s, method, margin, comp):
    if method in ('PLAPH', 'PLAPH-ind'):
        R = MAIN[s] if method == 'PLAPH' else IND[s]
        key = {('failure', 'linear'): ('beta_fp', 'beta_fn', 'mpl_beta_selected'),
               ('censoring', 'linear'): ('phi_fp', 'phi_fn', 'mpl_phi_selected'),
               ('failure', 'additive'): ('mpl_add_T_fp', 'mpl_add_T_fn', 'mpl_selected_additive_T'),
               ('censoring', 'additive'): ('mpl_add_C_fp', 'mpl_add_C_fn', 'mpl_selected_additive_C')}[(margin, comp)]
    elif method == 'PL':
        R = MAIN[s]
        key = {'linear': ('pl_fp', 'pl_fn', 'pl_beta_selected'),
               'additive': ('pl_add_fp', 'pl_add_fn', 'pl_selected_additive')}[comp]
    else:
        R = GRP[s]
        key = {'linear': ('lin_fp', 'lin_fn', 'lin_selected'),
               'additive': ('add_fp', 'add_fn', 'add_selected')}[comp]
    fp, fn, sz = key
    return [f3(exact(R, fp, fn).mean()), f3(R[fp].mean()), f3(R[fn].mean()), f'{R[sz].mean():.2f}']


# -------------------------------------------------------- linear coefficients
def coef_vals(s, method, margin, k):
    k -= 1
    if method in ('PLAPH', 'PLAPH-ind'):
        R = MAIN[s] if method == 'PLAPH' else IND[s]
        m = 'beta' if margin == 'failure' else 'phi'
        truth = (BETA if m == 'beta' else GAMMA)[k]
        est = R[f'mpl_{m}_active{k+1}'].to_numpy(float)
        ase = R[f'mpl_{m}_ase_active{k+1}'].to_numpy(float)
        cp = R[f'mpl_{m}_cp_active{k+1}'].to_numpy(float)
        cond = all4(R).to_numpy()
    elif method == 'PL':
        R = MAIN[s]; truth = BETA[k]
        est = R[f'pl_beta_active{k+1}'].to_numpy(float); ase = cp = None
        cond = (exact(R, 'pl_fp', 'pl_fn') & exact(R, 'pl_add_fp', 'pl_add_fn')).to_numpy()
    else:
        R = GRP[s]; truth = BETA[k]
        est = R[['beta1', 'beta3', 'beta6'][k]].to_numpy(float); ase = cp = None
        cond = (exact(R, 'lin_fp', 'lin_fn') & exact(R, 'add_fp', 'add_fn')).to_numpy()
    ok = np.isfinite(est)
    bias = est[ok].mean() - truth
    ese = est[ok].std(ddof=1)
    mse = np.mean((est[ok] - truth) ** 2)
    ec = est[cond & ok]
    ese_c = ec.std(ddof=1) if ec.size > 1 else np.nan
    if ase is None:
        return [f'{bias:.3f}', f'{ese:.3f}', '--', '--', f'{ese_c:.3f}', '--', '--', f'{mse:.4f}'], cond.sum()
    sel = ok & (est != 0) & np.isfinite(ase)
    ase_u = ase[sel].mean()
    cp_u = np.where(np.isfinite(cp), cp, 0.0)[ok].mean()
    selc = sel & cond
    ase_c = ase[selc].mean()
    cp_c = np.where(np.isfinite(cp), cp, 0.0)[ok & cond].mean()
    return [f'{bias:.3f}', f'{ese:.3f}', f'{ase_u:.3f}', f'{cp_u:.3f}',
            f'{ese_c:.3f}', f'{ase_c:.3f}', f'{cp_c:.3f}', f'{mse:.4f}'], cond.sum()


# --------------------------------------------------------------- the filler
src = open(tex).read()


def block(label):
    a = src.index(r'\label{' + label + '}')
    b = src.index(r'\end{', src.index(r'\bottomrule', a) if label.startswith('tab:supp-diag') else a + 10)
    return a, b


def fill(label, fn, end_token):
    global src
    a = src.index(r'\label{' + label + '}')
    b = src.index(end_token, a)
    lines = src[a:b].split('\n')
    n = 0
    for i, L in enumerate(lines):
        if r'\TBD{}' in L and re.match(r'^D\d+ ', L):
            new = fn(L)
            if new is not None:
                lines[i] = new; n += 1
    src = src[:a] + '\n'.join(lines) + src[b:]
    return n


def diag_fn(L):
    s = int(re.match(r'^D(\d+)', L).group(1))
    return diag_row(s)


def sel_fn(L):
    f = [x.strip() for x in L.rstrip(' \\').split('&')]
    s = int(f[0][1:]); method, margin, comp = f[1], f[2], f[3]
    if (method == 'grpsurv' and GRP[s] is None) or (method == 'PLAPH-ind' and IND[s] is None):
        return None
    return ' & '.join(f[:4] + sel_vals(s, method, margin, comp)) + r' \\'


cond_counts = []


def coef_fn_factory(margin):
    def fn(L):
        f = [x.strip() for x in L.rstrip(' \\').split('&')]
        s = int(f[0][1:]); method, k = f[1], int(f[2])
        if (method == 'grpsurv' and GRP[s] is None) or (method == 'PLAPH-ind' and IND[s] is None):
            return None
        vals, nc = coef_vals(s, method, margin, k)
        if method == 'PLAPH':
            cond_counts.append(nc)
        return ' & '.join(f[:4] + vals) + r' \\'
    return fn


print('diagnostics', fill('tab:supp-diagnostics', diag_fn, r'\end{tabular}'))
print('selection', fill('tab:e10-selection-detail', sel_fn, r'\end{longtable}'))
print('linear failure', fill('tab:e10-linear-failure', coef_fn_factory('failure'), r'\end{longtable}'))
print('linear censoring', fill('tab:e10-linear-censoring', coef_fn_factory('censoring'), r'\end{longtable}'))
print('conditional counts', min(cond_counts), max(cond_counts))
open(tex, 'w').write(src)

# pooled summaries for the main text (PLAPH, D1-D13)
out = []
for s in range(1, 14):
    for margin in ['failure', 'censoring']:
        for k in (1, 2, 3):
            v, _ = coef_vals(s, 'PLAPH', margin, k)
            out.append((s, margin, k, *[float(x) if x != '--' else np.nan for x in v]))
T = pd.DataFrame(out, columns=['s', 'margin', 'k', 'bias', 'ese', 'ase', 'cp', 'ese_c', 'ase_c', 'cp_c', 'mse'])
T['r'] = T.ase / T.ese; T['r_c'] = T.ase_c / T.ese_c
print(T.groupby('margin')[['cp', 'r', 'cp_c', 'r_c']].mean().round(3))
print(T.groupby(['s', 'margin']).cp_c.mean().unstack().round(3))
