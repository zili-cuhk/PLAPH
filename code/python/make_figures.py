"""Regenerate the article figures from the v2 outputs.

Usage: python3 make_figures.py SIM_RESULTS_DIR FIG_DIR [which ...]
(TCGA results are read from $PLAPH_TCGA_DIR, default SIM_RESULTS_DIR/tcga_new.)
which in {xi, luad, functions}; default all available.
"""
import glob, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RUNS, FIG = sys.argv[1], sys.argv[2]
WHICH = set(sys.argv[3:]) or {'xi', 'luad', 'functions'}
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
                     'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42})
CLIN = {'age', 'male', 'stage_II', 'stage_III', 'stage_IV'}


def save(fig, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, bbox_inches='tight'); plt.close(fig); print('wrote', path)


# ------------------------------------------------------------------- xi paths
def xi_paths():
    xis = [0, 5, 10, 15, 20, 25.6685, 30, 35, 40, 50]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True)
    for r, (base, n) in enumerate([(51, 500), (61, 700)]):
        frames = []
        for i, xi in enumerate(xis):
            d = glob.glob(f'{RUNS}/xi_sensitivity/V2_xi_sensitivity_S{base+i:02d}_*')[0]
            R = pd.concat([pd.read_csv(f) for f in glob.glob(d + '/chunk_*/scenario_001_replicates.csv')])
            frames.append(R[['seed', 'mpl_beta_active1', 'mpl_beta_inactive1']].assign(xi=xi))
        D = pd.concat(frames)
        for c, (col, truth, title) in enumerate([('mpl_beta_active1', 1.5, r'Signal $\beta_1=1.5$'),
                                                  ('mpl_beta_inactive1', 0.0, r'Null coordinate $\beta_2=0$')]):
            ax = axes[r, c]
            W = D.pivot(index='seed', columns='xi', values=col)
            for _, row in W.iterrows():
                ax.plot(xis, row.values, color=BLUE, alpha=0.06, lw=0.8)
            med = W.median()
            ax.plot(xis, med.values, color=BLUE, lw=2, marker='o', ms=4, label='Median')
            ax.axhline(truth, color=INK, ls='--', lw=1, label='Truth')
            ax.axvline(25.6685, color=MUTED, ls=':', lw=1)
            ax.set_title(f'{title}, $n={n}$', fontsize=9, color=INK)
            if r == 1: ax.set_xlabel(r'Kink coefficient $\xi$')
            if c == 0: ax.set_ylabel('Estimate')
            if r == 0 and c == 0: ax.legend(loc='upper right')
    fig.tight_layout()
    save(fig, f'{FIG}/supplement/E11_xi_coefficient_paths_n500_n700.pdf')


# ------------------------------------------------------------------- TCGA
def luad():
    tc = os.environ.get('PLAPH_TCGA_DIR', f'{RUNS}/tcga_new')
    fams = {'Frank': ([-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7], BLUE, 'o'),
            'Gaussian': ([-0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7], ORANGE, 's'),
            'Clayton': ([0.3, 0.5, 0.7], AQUA, '^'), 'Gumbel': ([0.3, 0.5, 0.7], YELLOW, 'D')}
    fig, axes = plt.subplots(3, 1, figsize=(5.6, 6.6), sharex=True)
    for fam, (taus, col, mk) in fams.items():
        ll, nT, nC = [], [], []
        for t in taus:
            gdir = f'{tc}/grid/{"Frank" if t == 0 else fam}_tau{t:+.2f}'
            c = pd.read_csv(gdir + '/ebic_choice.csv').iloc[0]
            ll.append(c.loglik); nT.append(c.rna_failure); nC.append(c.rna_censoring)
        jit = {'Frank': -0.012, 'Gaussian': 0.012, 'Clayton': -0.006, 'Gumbel': 0.006}[fam]
        x = np.array(taus) + jit
        axes[0].plot(x, ll, color=col, marker=mk, ms=5, lw=1.5, label=fam)
        axes[1].plot(x, nT, color=col, marker=mk, ms=5, lw=1.5)
        axes[2].plot(x, nC, color=col, marker=mk, ms=5, lw=1.5)
    axes[0].set_ylabel('Log likelihood'); axes[0].legend(ncol=4, loc='lower left', fontsize=8)
    axes[1].set_ylabel('Failure RNA\nselected'); axes[2].set_ylabel('Censoring RNA\nselected')
    for ax in axes[1:]: ax.set_yticks([0, 1, 2]); ax.set_ylim(-0.3, 2.3)
    axes[2].set_xlabel(r"Fixed Kendall's $\tau$ of the working copula")
    fig.tight_layout()
    save(fig, f'{FIG}/supplement/LUAD_fixed_tau_sensitivity.pdf')

    # subsampling frequencies
    def lin(d, m):
        L = pd.read_csv(f'{d}/primary_p500_q200_linear_{m}.csv')
        return set(L[(L.selected == 1) & ~L.feature.isin(CLIN)].feature)
    subs = sorted(glob.glob(f'{tc}/subsample/sub_*'))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    for ax, m, title in zip(axes, ['failure', 'censoring'], ['Failure margin', 'Censoring margin']):
        cnt = {}
        for d in subs:
            for g in lin(d, m): cnt[g] = cnt.get(g, 0) + 1
        top = sorted(cnt.items(), key=lambda x: -x[1])[:8][::-1]
        ax.barh([g for g, _ in top], [v / len(subs) for _, v in top], color=BLUE, height=0.6)
        for i, (_, v) in enumerate(top):
            ax.text(v / len(subs) + 0.01, i, f'{v/len(subs):.2f}', va='center', fontsize=8, color=MUTED)
        ax.axvline(0.5, color=MUTED, ls=':', lw=1)
        ax.set_xlim(0, 1); ax.set_xlabel('Selection frequency'); ax.set_title(title, fontsize=9)
        ax.grid(axis='y', visible=False)
    fig.tight_layout()
    save(fig, f'{FIG}/main/luad_subsampling_main.pdf')


# ---------------------------------------------------------------- functions
def function_figs():
    def curves(exp):
        f = f'{RUNS}/summary/curves_{exp}.csv'
        return pd.read_csv(f) if os.path.isfile(f) else None
    def grp_curves(s):
        fs = glob.glob(f'{RUNS}/grpsurv_results/V2_grpsurv_export_S{s:02d}_*/rows_*_curves.csv')
        if not fs: return None
        G = pd.concat([pd.read_csv(f) for f in fs])
        return G.groupby(['component', 'x']).fit.agg(mean='mean', lo=lambda v: v.quantile(.025),
                                                      hi=lambda v: v.quantile(.975)).reset_index()
    names = [r'$g_1(x)$', r'$g_2(x)$', r'$g_3(x)$', r'$g_4(x)$']
    M = curves('main')
    if M is None:
        print('curves_main.csv missing; skipped'); return
    sty = {'mpl_T': (BLUE, 'PLAPH'), 'pl': (ORANGE, 'PL'), 'grp': (AQUA, 'grpsurv'), 'ind': (YELLOW, 'PLAPH-ind')}

    def panel(ax, C, method, comp, color, label):
        c = C[(C.method == method) & (C.component == comp)].sort_values('x')
        ax.fill_between(c.x, c.lo, c.hi, color=color, alpha=0.15, lw=0)
        ax.plot(c.x, c['mean'], color=color, lw=1.6, label=label)
        return c

    # main text: failure margin, S4 and S6, PLAPH vs PL vs grpsurv
    for s in [4, 6]:
        fig, axes = plt.subplots(2, 2, figsize=(6.4, 4.8))
        G = grp_curves(s)
        for j, ax in enumerate(axes.flat, start=1):
            c = panel(ax, M[M.scenario_index == s], 'mpl_T', j, BLUE, 'PLAPH')
            panel(ax, M[M.scenario_index == s], 'pl', j, ORANGE, 'PL')
            if G is not None:
                g = G[G.component == j].sort_values('x')
                ax.plot(g.x, g['mean'], color=AQUA, lw=1.6, label='grpsurv')
            ax.plot(c.x, c.truth, color=INK, ls='--', lw=1, label='Truth')
            ax.set_title(names[j - 1], fontsize=9)
            if j == 1: ax.legend(fontsize=7)
        fig.tight_layout()
        save(fig, f'{FIG}/main/function_recovery_failure_S{s}.pdf')

    # supplement: combined failure + censoring, S1 and S9
    A = curves('ablation_tau0')
    for s in [1, 9]:
        fig, axes = plt.subplots(2, 4, figsize=(9.6, 4.6))
        G = grp_curves(s)
        for j in range(1, 5):
            ax = axes[0, j - 1]
            c = panel(ax, M[M.scenario_index == s], 'mpl_T', j, BLUE, 'PLAPH')
            panel(ax, M[M.scenario_index == s], 'pl', j, ORANGE, 'PL')
            if G is not None:
                g = G[G.component == j].sort_values('x'); ax.plot(g.x, g['mean'], color=AQUA, lw=1.6, label='grpsurv')
            if A is not None and (A.scenario_index == s).any():
                a = A[(A.scenario_index == s) & (A.method == 'mpl_T') & (A.component == j)].sort_values('x')
                ax.plot(a.x, a['mean'], color=YELLOW, lw=1.6, label='PLAPH-ind')
            ax.plot(c.x, c.truth, color=INK, ls='--', lw=1, label='Truth')
            ax.set_title('Failure ' + names[j - 1], fontsize=9)
            ax2 = axes[1, j - 1]
            c = panel(ax2, M[M.scenario_index == s], 'mpl_C', j, BLUE, 'PLAPH')
            if A is not None and (A.scenario_index == s).any():
                a = A[(A.scenario_index == s) & (A.method == 'mpl_C') & (A.component == j)].sort_values('x')
                ax2.plot(a.x, a['mean'], color=YELLOW, lw=1.6, label='PLAPH-ind')
            ax2.plot(c.x, c.truth, color=INK, ls='--', lw=1, label='Truth')
            ax2.set_title('Censoring ' + names[j - 1], fontsize=9)
        axes[0, 0].legend(fontsize=7); axes[1, 0].legend(fontsize=7)
        fig.tight_layout()
        save(fig, f'{FIG}/supplement/function_recovery_combined_S{s:02d}.pdf')

    # copula family F7 / F8
    F = curves('copula_family')
    if F is not None:
        for s, tag in [(37, 'F07'), (38, 'F08')]:
            fig, axes = plt.subplots(2, 4, figsize=(9.6, 4.6))
            for j in range(1, 5):
                for r, (meth, lab) in enumerate([('mpl_T', 'Failure '), ('mpl_C', 'Censoring ')]):
                    ax = axes[r, j - 1]
                    c = panel(ax, F[F.scenario_index == s], meth, j, BLUE, 'PLAPH')
                    ax.plot(c.x, c.truth, color=INK, ls='--', lw=1, label='Truth')
                    ax.set_title(lab + names[j - 1], fontsize=9)
            axes[0, 0].legend(fontsize=7)
            fig.tight_layout()
            save(fig, f'{FIG}/supplement/copula_family_function_recovery_{tag}.pdf')


if 'xi' in WHICH: xi_paths()
if 'luad' in WHICH: luad()
if 'functions' in WHICH: function_figs()
