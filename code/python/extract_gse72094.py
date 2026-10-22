"""Extract CYP17A1 expression and clinical fields from the GEO series matrix of GSE72094.

Download first (from NCBI GEO):
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72094/matrix/GSE72094_series_matrix.txt.gz
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15048&targ=self&form=text&view=data  (saved as GPL15048_table.txt)

Usage: python3 extract_gse72094.py GSE72094_series_matrix.txt.gz GPL15048_table.txt OUT.csv [GENE ...]
Writes one row per sample: clinical fields and every probe set mapped to the genes (default CYP17A1).
"""
import gzip, io, sys
import pandas as pd

matrix, platform, out = sys.argv[1:4]
genes = set(sys.argv[4:]) or {'CYP17A1'}
lines = open(platform).read().split('\n')
b = lines.index('!platform_table_begin')
e = lines.index('!platform_table_end') if '!platform_table_end' in lines else len(lines)
P = pd.read_csv(io.StringIO('\n'.join(lines[b + 1:e])), sep='\t', dtype=str)
probes = P[P.GeneSymbol.fillna('').str.split(' /// ').apply(lambda s: bool(genes & set(s)))]
print(probes[['ID', 'GeneSymbol', 'EntrezGeneID']].to_string(index=False))
ids = set(probes.ID)

keep = ('gender', 'age_at_diagnosis', 'vital_status', 'survival_time_in_days', 'Stage', 'patient_id', 'smoking_status')
ch, rows, header, intab = {}, [], None, False
with gzip.open(matrix, 'rt') as f:
    for line in f:
        if line.startswith('!Sample_geo_accession'):
            gsm = [x.strip('"') for x in line.rstrip('\n').split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            vals = [x.strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            keys = [v.split(': ')[0] for v in vals if v]
            if keys:
                k = max(set(keys), key=keys.count)
                if k in keep:
                    ch[k] = [v.split(': ', 1)[1] if ': ' in v else '' for v in vals]
        elif line.startswith('!series_matrix_table_begin'):
            intab = True; header = next(f)
        elif intab and not line.startswith('!'):
            if line.split('\t', 1)[0].strip('"') in ids:
                rows.append(line.rstrip('\n').split('\t'))
cols = [c.strip('"') for c in header.rstrip('\n').split('\t')]
E = pd.DataFrame(rows, columns=cols).set_index('ID_REF').T.astype(float)
E.columns = ['probe_' + c.strip('"') for c in E.columns]
D = pd.DataFrame(ch, index=gsm).join(E); D.index.name = 'gsm'
D.to_csv(out)
print(D.shape, 'written to', out)
