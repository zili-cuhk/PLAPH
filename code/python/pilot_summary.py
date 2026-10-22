"""Numbers of the 'Development pilot' paragraph (Supplement, tuning section).

Usage: python3 pilot_summary.py results/pilot
pilot_order4 / pilot_order3: cubic / quadratic basis on the earlier grid {1,1.5,2} x {0.10,...,0.30}
pilot_grid2: quadratic basis on the final grid {0.5,...,2.5} x {0.03,...,0.15}, budget 6000.
All on development seeds (root 26000000), designs D1, D4, D6, 20 replications each.
"""
import glob, sys
import warnings
import numpy as np, pandas as pd
warnings.simplefilter('ignore', pd.errors.PerformanceWarning)

root = sys.argv[1]
ERR = ['beta_fp', 'beta_fn', 'phi_fp', 'phi_fn', 'mpl_add_T_fp', 'mpl_add_T_fn', 'mpl_add_C_fp', 'mpl_add_C_fn']


def load(run, pat='replicates'):
    fs = glob.glob(f'{root}/{run}/spline_pilot/*/chunk_*/scenario_001_{pat}.csv')
    D = pd.concat([pd.read_csv(f).assign(scen=f.split('/')[-3][17:20]) for f in fs], ignore_index=True)
    return D.copy()


def exact_rate(R):
    pairs = [('beta_fp', 'beta_fn'), ('phi_fp', 'phi_fn'), ('mpl_add_T_fp', 'mpl_add_T_fn'), ('mpl_add_C_fp', 'mpl_add_C_fn')]
    return np.mean([((R[a] == 0) & (R[b] == 0)).mean() for a, b in pairs])


cub, quad, new = load('pilot_order4'), load('pilot_order3'), load('pilot_grid2')
print('basis comparison (earlier grid):')
print(f'  exact recovery averaged over the four supports: quadratic {exact_rate(quad):.3f}, cubic {exact_rate(cub):.3f}')
print(f'  acceptance: quadratic {quad.mpl_accepted.mean():.3f}, cubic {cub.mpl_accepted.mean():.3f}')
P = cub.merge(quad, on=['scen', 'seed'], suffixes=('_c', '_q'))
for m in ['T', 'C']:
    a, b = P[f'imse_mpl_add_{m}_c'], P[f'imse_mpl_add_{m}_q']
    print(f'  margin {m}: quadratic smaller in {np.mean(b < a):.0%}; median cubic/quadratic IMSE ratio {np.median(a / b):.2f}')
na = pd.concat([cub.assign(basis='cubic'), quad.assign(basis='quadratic')])
na = na[na.mpl_accepted != 1]
print(f'  non-accepted fits (both bases): {len(na)}')
print(na[['basis', 'seed', 'mpl_target_iterations', 'mpl_selection_gate_maximum', 'mpl_refit_stationarity_inf']].to_string(index=False))

print('earlier grid vs final grid (quadratic):')
for name, R in [('earlier', quad), ('final', new)]:
    ok = (R[ERR].sum(axis=1) == 0)
    print(f'  {name}: accepted {int(R.mpl_accepted.sum())}/60, all four supports exact {int(ok.sum())}/60, '
          f'max target iterations {int(R.mpl_target_iterations.max())}')
bad = quad[quad[ERR].sum(axis=1) > 0]
print('  earlier-grid misses:', bad[['scen', 'seed', 'beta_fn', 'phi_fn']].to_dict('records'))

T = load('pilot_grid2', 'tuning_path')
full, interior, only = 0, 0, []
for (s, sd), g in T.groupby(['scen', 'seed']):
    e = g[g.numerically_valid == 1]; m = e.criterion_value.min()
    tie = e[e.criterion_value <= m + max(1e-4, 1e-6 * abs(m))]
    full += len(tie) == len(e)
    inside = ((tie.c_L > 0.5) & (tie.c_L < 2.5) & (tie.kappa > 0.03) & (tie.kappa < 0.15)).any()
    interior += inside
    if not inside:
        only.append((s, sd, sorted(tie.c_L.unique()), sorted(tie.kappa.unique())))
print(f'final grid: selected model optimal at all grid points in {full}/60; an interior point optimal in {interior}/60')
print('  boundary-only replications:', only)
