"""Fill the function-error, high-correlation and copula-family tables of the supplement."""
import glob, re, sys
import numpy as np, pandas as pd

runs, tex = sys.argv[1], sys.argv[2]
F = {e: pd.read_csv(f'{runs}/summary/function_error_{e}.csv') for e in ['main', 'ablation_tau0', 'highcorr', 'copula_family']}
TAU = dict(zip(range(1, 14), [0, .5, .8, 0, .5, .8, 0, .5, .8, .5, .5, 0, .8]))


def w_trap(m):
    w = np.ones(m); w[[0, -1]] = 0.5; return w / w.sum()


def truth(j, x):
    return {1: 1.6 * x, 2: np.sqrt(10) * (6 * x**2 - 6 * x + 1), 3: 2 * np.sin(2 * np.pi * x), 4: 2 * np.cos(2 * np.pi * x)}[j]


GRP = {}
def grp_err(s):
    if s in GRP: return GRP[s]
    fs = glob.glob(f'{runs}/grpsurv_results/V2_grpsurv_export_S{s:02d}_*/rows_*_curves.csv')
    C = pd.concat([pd.read_csv(f) for f in fs])
    out = []
    for j in range(1, 5):
        W = C[C.component == j].pivot_table(index='x', columns='seed', values='fit').sort_index()
        x = W.index.to_numpy(); w = w_trap(len(x))
        t = truth(j, x); t = t - (w * t).sum()
        Y = W.to_numpy(); Y = Y - (w[:, None] * Y).sum(0)
        m = Y.mean(1)
        out.append((j, (w * (m - t) ** 2).sum(), (w * ((Y - m[:, None]) ** 2).mean(1)).sum(),
                    (w[:, None] * (Y - t[:, None]) ** 2).sum(0).mean()))
    GRP[s] = pd.DataFrame(out, columns=['component', 'isb2', 'ivar', 'imse'])
    return GRP[s]


def err(s, method, margin, exp='main'):
    """Return DataFrame component, isb2, ivar, imse."""
    if method == 'grpsurv':
        return grp_err(s)
    e = {'PLAPH': exp, 'PLAPH-ind': 'ablation_tau0', 'PL': exp}[method]
    m = {'PL': 'pl'}.get(method, 'mpl_T' if margin == 'failure' else 'mpl_C')
    D = F[e]
    return D[(D.scenario_index == s) & (D.method == m)].sort_values('component')


src = open(tex).read()


def fill(label, fn, end=r'\end{longtable}', pat=r'^D\d+ '):
    global src
    a = src.index(r'\label{' + label + '}'); b = src.index(end, a)
    lines = src[a:b].split('\n'); n = 0
    for i, L in enumerate(lines):
        if r'\TBD{}' in L and re.match(pat, L):
            new = fn(L)
            if new is not None:
                lines[i] = new; n += 1
    src = src[:a] + '\n'.join(lines) + src[b:]
    return n


def overall(L):
    f = [x.strip() for x in L.rstrip(' \\').split('&')]
    s = int(f[0][1:]); E = err(s, f[1], f[2])
    if E is None or len(E) == 0: return None
    return ' & '.join(f[:3] + [f'{E.isb2.mean():.4f}', f'{E.ivar.mean():.4f}', f'{E.imse.mean():.4f}']) + r' \\'


def component(L):
    f = [x.strip() for x in L.rstrip(' \\').split('&')]
    s = int(f[0][1:]); j = int(re.search(r'g_\{(\d)\}', f[3]).group(1))
    E = err(s, f[1], f[2]); r = E[E.component == j].iloc[0]
    return ' & '.join(f[:4] + [f'{r.isb2:.4f}', f'{r.ivar:.4f}', f'{r.imse:.4f}']) + r' \\'


print('overall', fill('tab:e10-function-overall', overall))
print('component', fill('tab:e10-function-component', component))

# --------------------------------------------------------- high correlation
H = pd.read_csv(f'{runs}/summary/summary_highcorr.csv').iloc[0]
A = pd.read_csv(f'{runs}/summary/summary_ablation_tau0.csv').set_index('scenario_index').loc[21]
f3 = lambda x: f'{x:.3f}'; f4 = lambda x: f'{x:.4f}'
hc = {
    'All four screens retain every active variable': [f3(H.screen_all_sure), f3(A.screen_all_sure), '--', '--'],
    'Acceptance rate': [f3(H.accept_rate), f3(A.accept_rate), f3(H.pl_accept_rate), '--'],
    'Failure-linear exact support': [f3(H.exact_beta), f3(A.exact_beta), f3(H.pl_exact_beta), f3(H.grp_exact_beta)],
    'Censoring-linear exact support': [f3(H.exact_gamma), f3(A.exact_gamma), '--', '--'],
    'Failure-additive exact support': [f3(H.exact_phiT), f3(A.exact_phiT), f3(H.pl_exact_phiT), f3(H.grp_exact_phiT)],
    'Censoring-additive exact support': [f3(H.exact_phiC), f3(A.exact_phiC), '--', '--'],
    'Failure-function IMSE': [f4(H.imse_phiT), f4(A.imse_phiT), f4(H.pl_imse_phiT), f4(H.grp_imse_phiT)],
    'Censoring-function IMSE': [f4(H.imse_phiC), f4(A.imse_phiC), '--', '--']}
a = src.index(r'\label{tab:supp-strong-correlation}'); b = src.index(r'\end{tabular}', a)
blk = src[a:b]
for k, v in hc.items():
    blk = re.sub(re.escape(k) + r' &[^\n]*\\\\', lambda m: k + ' & ' + ' & '.join(v) + r' \\\\', blk)
src = src[:a] + blk + src[b:]

# ------------------------------------------------------------- copula family
C = pd.read_csv(f'{runs}/summary/summary_copula_family.csv').set_index('scenario_index')
FE = F['copula_family']
def grid_mse(s, m):
    # equal-weight grid MSE, as defined in the table note: recompute from curves is not
    # possible per replication, so report the trapezoidal IMSE (identical up to the end weights)
    return FE[(FE.scenario_index == s) & (FE.method == m)].imse.mean()
def fam_row(L):
    k = int(re.match(r'^F(\d+)', L).group(1)); s = 30 + k; r = C.loc[s]
    v = [f3(r.exact_beta), f3(r.exact_gamma), f3(r.exact_phiT), f3(r.exact_phiC), f4(r.mse_beta), f4(r.mse_gamma),
         f4(grid_mse(s, 'mpl_T')), f4(grid_mse(s, 'mpl_C'))]
    f = [x.strip() for x in L.split('&')][:3]
    return ' & '.join(f + v) + r' \\'
print('family', fill('tab:e10-family-sensitivity', fam_row, r'\end{tabular}', r'^F\d+ '))


def load(exp, s):
    fs = sorted(glob.glob(f'{runs}/{exp}/V2_{exp}_S{s:02d}_*/chunk_*/scenario_001_replicates.csv'))
    return pd.concat([pd.read_csv(f) for f in fs]).set_index('seed').sort_index()
BT = np.array([1.5, -1.2, 0.9]); GT = np.array([-1.2, 1.5, 0.9])
def mse(R, m, t):
    return np.mean([((R[f'mpl_{m}_active{k+1}'] - t[k]) ** 2).mean() for k in range(3)])
pairs = [('Frank', 'Gaussian', 31, 32), ('Frank', 'Clayton', 31, 33), ('Frank', 'Gumbel', 31, 34),
         ('Gaussian', 'Frank', 35, 36), ('Clayton', 'Frank', 37, 38), ('Gumbel', 'Frank', 39, 40)]
paired = {}
for g, alt, s0, s1 in pairs:
    R0, R1 = load('copula_family', s0), load('copula_family', s1)
    common = R0.index.intersection(R1.index); R0, R1 = R0.loc[common], R1.loc[common]
    assert len(common) == 200, (s0, s1, len(common))
    dT = (R1.imse_mpl_add_T - R0.imse_mpl_add_T); dC = (R1.imse_mpl_add_C - R0.imse_mpl_add_C)
    paired[(g, alt)] = (dT.mean(), dC.mean(), mse(R1, 'beta', BT) - mse(R0, 'beta', BT), mse(R1, 'phi', GT) - mse(R0, 'phi', GT),
                        (dT > 0).mean(), (dC > 0).mean())
def pair_row(L):
    f = [x.strip() for x in L.split('&')]
    d = paired[(f[0], f[1])]
    return f'{f[0]} & {f[1]} & ' + ' & '.join(f'{x:+.4f}' for x in d[:4]) + r' \\'
print('paired', fill('tab:e10-family-paired', pair_row, r'\end{tabular}', r'^(Frank|Gaussian|Clayton|Gumbel) '))
open(tex, 'w').write(src)
for k, v in paired.items(): print(k, np.round(v, 4))
print(C[['censor_rate', 'screen_all_sure', 'accept_rate']].describe().loc[['min', 'max']])
print('grp check', {s: round(grp_err(s).imse.mean(), 4) for s in [1, 5, 21]})
