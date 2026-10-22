"""Fill the remaining supplement tables whose bodies are generated wholesale:
kink-coefficient path, TCGA gamma = 0.5 re-selection, TCGA additive-support
sensitivity, TCGA screened-dimension sensitivity and the external evaluation.

Usage: python3 fill_supp_extra.py SIM_RESULTS_DIR TCGA_RESULTS_DIR EXTERNAL_JSON SUPPLEMENT_TEX
"""
import glob, json, os, re, sys
import numpy as np, pandas as pd

runs, tc, ext_json, tex = sys.argv[1:5]
CLIN = {'age', 'male', 'stage_II', 'stage_III', 'stage_IV'}
src = open(tex).read()


def body(label, rows):
    """Replace everything between \\midrule and \\bottomrule of the table with this label."""
    global src
    a = src.index(r'\label{' + label + '}')
    m = src.index(r'\midrule', a) + len(r'\midrule'); e = src.index(r'\bottomrule', a)
    src = src[:m] + '\n' + '\n'.join(rows) + '\n' + src[e:]


def neg(x, fmt):
    v = format(x, fmt)
    return f'$-${v[1:]}' if v.startswith('-') else v


# ------------------------------------------------------------ kink coefficient
xis = [0, 5, 10, 15, 20, 25.6685, 30, 35, 40, 50]
def xi_stats(s):
    d = glob.glob(f'{runs}/xi_sensitivity/V2_xi_sensitivity_S{s:02d}_*')[0]
    x = pd.concat([pd.read_csv(f) for f in glob.glob(d + '/chunk_*/scenario_001_replicates.csv')])
    j = ((x.beta_fp == 0) & (x.beta_fn == 0) & (x.phi_fp == 0) & (x.phi_fn == 0)).mean()
    return f"{x.mpl_beta_active1.median():.3f} & {(x.mpl_beta_inactive1 == 0).mean():.3f} & {j:.3f}"
body('tab:e11-xi-path', [f"{x:g} & {xi_stats(51 + i)} & {xi_stats(61 + i)} \\\\" for i, x in enumerate(xis)])


# ------------------------------------------------------- TCGA helpers
def lin(d, m):
    L = pd.read_csv(f'{d}/primary_p500_q200_linear_{m}.csv')
    L = L[(L.selected == 1) & ~L.feature.isin(CLIN)]
    return dict(zip(L.feature, L.coefficient))


def choice(fam, tau, gamma):
    g = f'{tc}/grid/{fam}_tau{tau:+.2f}'
    T = pd.read_csv(f'{g}/ebic_grid.csv'); T['ms'] = (T.ebic - T.bic) / 2
    v = T[T.accepted == 1].copy(); v['crit'] = v.bic + 2 * gamma * v.ms
    m = v.crit.min(); tie = v[v.crit <= m + max(1e-4, 1e-6 * abs(m))]
    b = tie.sort_values(['c_L', 'kappa'], ascending=False).iloc[0]
    return b, f'{g}/{os.path.basename(b.output_dir)}'


# gamma = 0.5 re-selection (Frank and Gaussian rows)
rows = []
for fam, taus in [('Frank', [-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7]), ('Gaussian', [-0.7, -0.5, -0.3, 0.3, 0.5, 0.7])]:
    if fam == 'Gaussian': rows.append(r'\addlinespace')
    for t in taus:
        b, d = choice(fam, t, 0.5); T, C = lin(d, 'failure'), lin(d, 'censoring')
        genes = ', '.join(f'{g} ({neg(c, ".2f")})' for g, c in sorted(T.items(), key=lambda z: z[1]))
        name = 'Independence' if t == 0 else fam
        rows.append(f"{name} & ${t:g}$ & {b.c_L:g} & {neg(b.loglik, '.1f')} & {len(T)}/{len(C)} & {genes} \\\\")
body('tab:luad-gamma-half', rows)


# additive-support sensitivity
def add_row(label, d, kappa):
    r = pd.read_csv(d + '/fit_row.csv').iloc[0]
    return (f"{label} & {kappa:g} & {kappa:g} & {int(r.rna_failure)}/{int(r.rna_censoring)} & "
            f"{int(r.cna_failure)}/{int(r.cna_censoring)} & {int(r.df)} & {neg(r.loglik, '.1f')} & "
            f"{r.bic:.1f} & {r.ebic:.1f} & {'yes' if r.accepted == 1 else 'no'} \\\\")
rows = [add_row('EBIC-selected fit', f'{tc}/grid/Frank_tau+0.00/cL6.00_kappa0.15', 0.15)]
rows += [add_row(f'$\\kappa={k:g}$', f'{tc}/additive/kpath_{k:.4f}', k) for k in [0.02, 0.015, 0.01, 0.0075, 0.005, 0.0025]]
rows.append(add_row('Weight exponent 0.25', f'{tc}/additive/weight_exponent_0.25', 0.15))
body('tab:luad-corrected-additive-audit', rows)


# screened-dimension sensitivity
def dim_row(pq, d):
    c = pd.read_csv(d + '/ebic_choice.csv').iloc[0]
    return (f"$({pq[0]},{pq[1]})$ & {int(c.rna_failure)} & {int(c.rna_censoring)} & {int(c.cna_failure)} & "
            f"{int(c.cna_censoring)} & {neg(c.loglik, '.1f')} & {c.bic:.1f} & {c.ebic:.1f} \\\\")
body('tab:luad-dimension-audit', [dim_row((500, 200), f'{tc}/grid/Frank_tau+0.00'),
                                  dim_row((250, 200), f'{tc}/dimension/p250_q200'),
                                  dim_row((500, 100), f'{tc}/dimension/p500_q100')])


# external evaluation
E = json.load(open(ext_json))
rows = [f"{g} & {neg(v['tcga'], '.3f')} & {neg(v['external'], '.3f')} ({v['se']:.3f}) & {'yes' if v['sign_agree'] else 'no'} & {v['p']:.2f} \\\\"
        for g, v in E['genes'].items()]
body('tab:external-supp', rows)
open(tex, 'w').write(src)
print('extra tables written')
