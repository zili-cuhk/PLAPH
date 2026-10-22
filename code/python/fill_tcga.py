"""Fill TCGA-LUAD tables/text (gamma = 1 primary; gamma = 0.5 sensitivity)."""
import glob, os, sys, re
import pandas as pd, numpy as np
tcga, main_tex, supp_tex = sys.argv[1:4]
CLIN = {'age', 'male', 'stage_II', 'stage_III', 'stage_IV'}

def lin(d, m):
    L = pd.read_csv(f'{d}/primary_p500_q200_linear_{m}.csv')
    L = L[(L.selected == 1) & ~L.feature.isin(CLIN)]
    return dict(zip(L.feature, L.coefficient))

def add(d, m):
    A = pd.read_csv(f'{d}/primary_p500_q200_additive_{m}.csv')
    return set(A[A.selected == 1].feature)

def grid_dir(fam, tau):
    return f'{tcga}/grid/{fam}_tau{tau:+.2f}'

def choice(fam, tau, gamma=1.0):
    T = pd.read_csv(f'{grid_dir(fam, tau)}/ebic_grid.csv')
    T['ms'] = (T.ebic - T.bic) / 2
    v = T[T.accepted == 1].copy(); v['crit'] = v.bic + 2 * gamma * v.ms
    m = v.crit.min(); tie = v[v.crit <= m + max(1e-4, 1e-6 * abs(m))]
    b = tie.sort_values(['c_L', 'kappa'], ascending=False).iloc[0]
    d = f'{grid_dir(fam, tau)}/{os.path.basename(b.output_dir)}'
    return b, d

N, D = 482, 174; NC = N - D
lamT = lambda c: c * np.log(D) / N
lamC = lambda c: c * np.log(NC) / N
ref_d = choice('Frank', 0)[1]
refT, refC = set(lin(ref_d, 'failure')), set(lin(ref_d, 'censoring'))
jac = lambda a, b: 1.0 if not a and not b else len(a & b) / len(a | b)
fmt_j = lambda x: f'{x:.2f}' if x not in (0, 1) else str(int(x))

main = open(main_tex).read(); supp = open(supp_tex).read()

def fill_rows(src, label, fn, end=r'\end{tabular}'):
    a = src.index(r'\label{' + label + '}'); b = src.index(end, a)
    lines = src[a:b].split('\n'); n = 0
    for i, L in enumerate(lines):
        if r'\TBD{}' in L and re.match(r'^(Frank|Gaussian|Independence|Clayton|Gumbel) ', L):
            lines[i] = fn(L); n += 1
    return src[:a] + '\n'.join(lines) + src[b:], n

def parse(L):
    f = [x.strip() for x in L.split('&')]
    fam = f[0]; tau = float(f[1].strip('$'))
    return ('Frank' if fam == 'Independence' else fam), tau, f

def main_row(L):
    fam, tau, f = parse(L)
    b, d = choice(fam, tau)
    T, C = lin(d, 'failure'), lin(d, 'censoring')
    aT, aC = add(d, 'failure'), add(d, 'censoring')
    jt = '1' if tau == 0 else fmt_j(jac(set(T), refT)); jc = '1' if tau == 0 else fmt_j(jac(set(C), refC))
    return f'{f[0]} & {f[1]} & {b.loglik:.1f} & {len(T)}/{len(C)} & {len(aT)}/{len(aC)} & {jt} & {jc} \\\\'

def supp_row(L):
    fam, tau, f = parse(L)
    b, d = choice(fam, tau)
    return (f'{f[0]} & {f[1]} & $({b.c_L:g},{b.kappa:g})$ & {lamT(b.c_L):.4f} & {lamC(b.c_L):.4f} & '
            f'{b.ebic:.1f} & ${b.refit_residual:.1e}'.replace('e-0', r'\times10^{-').replace('e-', r'\times10^{-') + '}$ \\\\')

main, n1 = fill_rows(main, 'tab:luad-grid-main', main_row)
supp, n2 = fill_rows(supp, 'tab:luad-fixed-supp', supp_row)
print('rows', n1, n2)

# support table over the 13 Frank/Gaussian fits
fits = [('Frank', t) for t in [-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7]] + [('Gaussian', t) for t in [-0.7, -0.5, -0.3, 0.3, 0.5, 0.7]]
coef = {'failure': {}, 'censoring': {}}
for fam, t in fits:
    d = choice(fam, t)[1]
    for m in coef:
        for g, c in lin(d, m).items():
            coef[m].setdefault(g, []).append(c)
ref = {'failure': lin(ref_d, 'failure'), 'censoring': lin(ref_d, 'censoring')}
rows = []
for m in ['failure', 'censoring']:
    for g, cs in sorted(coef[m].items(), key=lambda x: -len(x[1])):
        c0 = f'{ref[m][g]:.3f}' if g in ref[m] else '--'
        rows.append(f'{m.capitalize()} & {g} & {c0} & {len(cs)} & ${min(cs):.3f}$ to ${max(cs):.3f}$ \\\\')
a = supp.index(r'\label{tab:luad-support-full}'); s0 = supp.index(r'\midrule', a) + len(r'\midrule'); e0 = supp.index(r'\bottomrule', a)
supp = supp[:s0] + '\n' + '\n'.join(rows) + '\n' + supp[e0:]

# profile table
P = {r.family: r for r in (pd.read_csv(f).iloc[0] for f in glob.glob(f'{tcga}/profile/*/fit_row.csv'))}
def prof_row(L):
    fam = L.split('&')[0].strip(); r = P[fam]
    return (f'{fam} & {r.tau_hat:.3f} & {r.loglik:.1f} & {int(r.df)} & {r.aic:.1f} & {r.bic:.1f} & '
            f'{int(r.rna_failure)}/{int(r.rna_censoring)} & {int(r.cna_failure)}/{int(r.cna_censoring)} \\\\')
supp, n3 = fill_rows(supp, 'tab:luad-profile', prof_row)
print('profile rows', n3)

# subsamples (Frank tau = 0, full-data tuning)
fT, fC, sizes = {}, {}, []
for d in sorted(glob.glob(f'{tcga}/subsample/sub_*')):
    T, C = set(lin(d, 'failure')), set(lin(d, 'censoring'))
    for g in T: fT[g] = fT.get(g, 0) + 1
    for g in C: fC[g] = fC.get(g, 0) + 1
    r = pd.read_csv(f'{d}/fit_row.csv').iloc[0]
    sizes.append((len(T), len(C), len(add(d, 'failure')), len(add(d, 'censoring')), jac(T, refT), jac(C, refC), r.accepted, r.n, r.deaths))
S = pd.DataFrame(sizes, columns=['nT', 'nC', 'aT', 'aC', 'jT', 'jC', 'acc', 'n', 'deaths'])
R = len(S)
rows = [f'Failure & {g} & {fT.get(g, 0)/R:.3f} \\\\' for g in sorted(refT)] + [r'\midrule'] + \
       [f'Censoring & {g} & {fC.get(g, 0)/R:.3f} \\\\' for g in sorted(refC)]
a = supp.index(r'\label{tab:luad-subsampling-full-support}'); s0 = supp.index(r'\midrule', a) + len(r'\midrule'); e0 = supp.index(r'\bottomrule', a)
supp = supp[:s0] + '\n' + '\n'.join(rows) + '\n' + supp[e0:]
open(main_tex, 'w').write(main); open(supp_tex, 'w').write(supp)
print(S.mean().round(3).to_dict())
print('top T', sorted(fT.items(), key=lambda x: -x[1])[:5], 'top C', sorted(fC.items(), key=lambda x: -x[1])[:5])

# gamma = 0.5 sensitivity summary (printed for the text and the supplementary table)
for fam, t in fits:
    b, d = choice(fam, t, 0.5)
    print('g0.5', fam, t, b.c_L, sorted(lin(d, 'failure')), len(lin(d, 'censoring')))
