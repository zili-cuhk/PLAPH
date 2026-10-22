"""Regenerate every numeric table of the article from the result files and
compare it with the values reported in the article.

Usage:
  python3 reproduce_tables.py SIM_RESULTS_DIR TCGA_RESULTS_DIR EXTERNAL_JSON MAIN_TABLES SUPP_TABLES OUT_DIR

MAIN_TABLES and SUPP_TABLES are reference_tables/main_tables.tex and
supp_tables.tex: the numeric tables of the main article and the Supplementary
Material as LaTeX (captions and notes removed). They are only read. Their rows
supply the row keys (setting, method, margin, ...). All value cells are blanked
and refilled by the fill_*.py scripts in this directory from the result files;
the regenerated table bodies are written to OUT_DIR/tables/<label>.tex and every
cell whose number differs from the article is listed in OUT_DIR/table_check.txt.
"""
import os, re, shutil, subprocess, sys

res, tcga, ext_json, main_tex, supp_tex, out = sys.argv[1:7]
here = os.path.dirname(os.path.abspath(__file__))
work = os.path.join(out, 'work'); os.makedirs(work, exist_ok=True)
os.makedirs(os.path.join(out, 'tables'), exist_ok=True)
M, S = os.path.join(work, 'main.tex'), os.path.join(work, 'supp.tex')
shutil.copy(main_tex, M); shutil.copy(supp_tex, S)

# label -> (file, number of key columns, end token); None key count = body regenerated wholesale
TABLES = {
    'tab:e10-selection': ('main', 1, r'\end{tabular}'), 'tab:e10-performance': ('main', 1, r'\end{tabular}'),
    'tab:dimension-main': ('main', 1, r'\end{tabular}'), 'tab:luad-grid-main': ('main', 2, r'\end{tabular}'),
    'tab:supp-diagnostics': ('supp', 1, r'\end{tabular}'), 'tab:e10-selection-detail': ('supp', 4, r'\end{longtable}'),
    'tab:e10-linear-failure': ('supp', 4, r'\end{longtable}'), 'tab:e10-linear-censoring': ('supp', 4, r'\end{longtable}'),
    'tab:e10-function-overall': ('supp', 3, r'\end{longtable}'), 'tab:e10-function-component': ('supp', 4, r'\end{longtable}'),
    'tab:supp-strong-correlation': ('supp', 1, r'\end{tabular}'), 'tab:e10-family-sensitivity': ('supp', 3, r'\end{tabular}'),
    'tab:e10-family-paired': ('supp', 2, r'\end{tabular}'), 'tab:luad-fixed-supp': ('supp', 2, r'\end{tabular}'),
    'tab:luad-profile': ('supp', 1, r'\end{tabular}'),
    'tab:luad-support-full': ('supp', None, r'\end{tabular}'), 'tab:luad-subsampling-full-support': ('supp', None, r'\end{tabular}'),
    'tab:luad-gamma-half': ('supp', None, r'\end{tabular}'), 'tab:e11-xi-path': ('supp', None, r'\end{tabular}'),
    'tab:luad-corrected-additive-audit': ('supp', None, r'\end{tabular}'), 'tab:luad-dimension-audit': ('supp', None, r'\end{tabular}'),
    'tab:external-supp': ('supp', None, r'\end{tabular}'),
}


def rows_of(text, label, end):
    a = text.index(r'\label{' + label + '}')
    m = text.index(r'\midrule', a); e = text.index(end, a)
    if 'longtable' in end:                 # rows follow \endlastfoot in a longtable
        return a, m, e
    b = text.rfind(r'\bottomrule', a, e)
    return a, m, (b if b > 0 else e)


def is_row(L):
    L = L.strip()
    return '&' in L and L.endswith('\\\\') and not L.startswith(('Setting', 'Family', 'Gene', 'Margin', 'Configuration', r'\multicolumn', '$\\tau$', 'Metric', 'Generating', r'\(\xi', '&'))


def blank(path):
    t = open(path).read()
    for lab, (f, k, end) in TABLES.items():
        if k is None or ('main' if path == M else 'supp') != f: continue
        a, m, b = rows_of(t, lab, end)
        lines = t[m:b].split('\n')
        for i, L in enumerate(lines):
            if is_row(L) and not L.lstrip().startswith(r'\endhead'):
                cells = L.rstrip().rstrip('\\').split('&')
                lines[i] = '&'.join(cells[:k] + [r' \TBD{} '] * (len(cells) - k)) + r'\\'
        t = t[:m] + '\n'.join(lines) + t[b:]
    open(path, 'w').write(t)


blank(M); blank(S)
py = sys.executable
run = lambda *a: subprocess.run([py, os.path.join(here, a[0]), *a[1:]], check=True, capture_output=True, text=True)
run('fill_main_sim.py', res, M)
run('fill_supp_sim.py', res, S)
run('fill_supp_functions.py', res, S)
run('fill_tcga.py', tcga, M, S)
run('fill_supp_extra.py', res, tcga, ext_json, S)

NUM = re.compile(r'-?\d+(?:\.\d+)?(?:\\times10\^\{-?\d+\})?')
def nums(s):
    s = s.replace('$-$', '-').replace('$', '').replace(r'\(', '').replace(r'\)', '')
    return NUM.findall(s)

orig = {'main': open(main_tex).read(), 'supp': open(supp_tex).read()}
new = {'main': open(M).read(), 'supp': open(S).read()}
report, n_cells = [], 0
for lab, (f, k, end) in TABLES.items():
    _, m0, b0 = rows_of(orig[f], lab, end); _, m1, b1 = rows_of(new[f], lab, end)
    r0 = [L for L in orig[f][m0:b0].split('\n') if is_row(L)]
    r1 = [L for L in new[f][m1:b1].split('\n') if is_row(L)]
    open(os.path.join(out, 'tables', lab.replace(':', '_') + '.tex'), 'w').write('\n'.join(r1) + '\n')
    if len(r0) != len(r1):
        report.append(f'{lab}: {len(r0)} reference rows vs {len(r1)} regenerated rows'); continue
    for L0, L1 in zip(r0, r1):
        a, b = nums(L0), nums(L1); n_cells += len(a)
        if a != b:
            report.append(f'{lab}:\n  article    : {L0.strip()}\n  regenerated: {L1.strip()}')
    if r'\TBD' in new[f][m1:b1]:
        report.append(f'{lab}: unfilled cells remain')
with open(os.path.join(out, 'table_check.txt'), 'w') as fh:
    fh.write(f'{len(TABLES)} tables, {n_cells} numeric cells compared, {len(report)} discrepancies\n\n' + '\n'.join(report) + '\n')
print(f'{len(TABLES)} tables, {n_cells} numeric cells compared, {len(report)} discrepancies (see {out}/table_check.txt)')
