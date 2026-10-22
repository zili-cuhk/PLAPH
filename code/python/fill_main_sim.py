"""Fill the simulation tables of PLAPH_main.tex from the v2 summaries."""
import pandas as pd, re, sys
runs, tex = sys.argv[1], sys.argv[2]
S = pd.read_csv(f'{runs}/summary/summary_main.csv').set_index('scenario_index')
A = pd.read_csv(f'{runs}/summary/summary_ablation_tau0.csv').set_index('scenario_index')
P = pd.read_csv(f'{runs}/summary/ablation_paired.csv')
pr = P[P.metric == 'imse_mpl_add_T'].set_index('scenario_index').share_favoring_true_tau
f3 = lambda x: f'{x:.3f}'
f4 = lambda x: f'{x:.4f}'

def t1(s):
    r = S.loc[s]
    ind = ['--'] * 4 if r.tau == 0 else [f3(A.loc[s, c]) for c in ['exact_beta', 'exact_gamma', 'exact_phiT', 'exact_phiC']]
    v = [f3(r[c]) for c in ['exact_beta', 'exact_gamma', 'exact_phiT', 'exact_phiC']] + ind + \
        [f3(r.pl_exact_beta), f3(r.pl_exact_phiT), f3(r.grp_exact_beta), f3(r.grp_exact_phiT)]
    return f'{r.tau:.2f} & ' + ' & '.join(v) + r' \\'

def t2(s):
    r = S.loc[s]; z = r.tau == 0
    g = lambda c: '--' if z else f4(A.loc[s, c])
    v = [f4(r.mse_beta), g('mse_beta'), f4(r.pl_mse_beta), f4(r.grp_mse_beta),
         f4(r.mse_gamma), g('mse_gamma'),
         f4(r.imse_phiT), g('imse_phiT'), f4(r.pl_imse_phiT), f4(r.grp_imse_phiT),
         f4(r.imse_phiC), g('imse_phiC'), '--' if z else f3(pr.loc[s])]
    return f'{r.tau:.2f} & ' + ' & '.join(v) + r' \\'

def t3(s):
    r = S.loc[s]
    v = [f3(r[c]) for c in ['exact_beta', 'exact_gamma', 'exact_phiT', 'exact_phiC']] + \
        [f4(r.imse_phiT), f4(A.loc[s, 'imse_phiT']), f4(r.pl_imse_phiT), f4(r.grp_imse_phiT)]
    return f'{int(r.p)} & ' + ' & '.join(v) + r'\\'

src = open(tex).read()
def fill_block(src, label, rows):
    a = src.index(r'\label{' + label + '}')
    b = src.index(r'\end{tabular}', a)
    blk = src[a:b]
    lines = blk.split('\n')
    k = 0
    for i, L in enumerate(lines):
        if r'\TBD{}' in L and re.match(r'^\s*[0-9]', L):
            lines[i] = rows[k]; k += 1
    assert k == len(rows), (label, k, len(rows))
    return src[:a] + '\n'.join(lines) + src[b:]
order = [1, 2, 3, 4, 5, 6, 7, 8, 9]
src = fill_block(src, 'tab:e10-selection', [t1(s) for s in order])
src = fill_block(src, 'tab:e10-performance', [t2(s) for s in order])
src = fill_block(src, 'tab:dimension-main', [t3(s) for s in [10, 5, 11]])
open(tex, 'w').write(src)
print('filled')
