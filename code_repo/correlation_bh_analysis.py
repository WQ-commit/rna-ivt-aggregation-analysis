# -*- coding: utf-8 -*-
"""
Correlation and Benjamini-Hochberg (BH) analysis for Fig. 1c and Fig. 1f.

Inputs (data/):
  - Fig1f_plot_data_n9.csv : 9 diverse transcripts, 8 sequence/structure
    parameters x 4 aggregation metrics (Fig. 1f).
  - transcripts_n13.csv    : 13 transcripts with onset concentration vs length (Fig. 1c).

Outputs (output/):
  - Fig1f_r_matrix.csv, Fig1f_p_matrix.csv, Fig1f_q_matrix.csv

Run:  python correlation_bh_analysis.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
OUT = os.path.join(HERE, 'output')
os.makedirs(OUT, exist_ok=True)

PARAMS = ['Length_nt', 'GC', 'MFE_per_nt', 'Pair_frac', 'Stems', 'Max_loop', 'LZ76', 'k2_H']
METRICS = ['Onset_time_min', 'Onset_conc_ng_uL', 'Rate_mass_ng_uL_min', 'Rate_molar_pM_min']
LABELS = ['Onset time', 'Onset conc.', 'Rate (mass)', 'Rate (molar)']


def bh(pvals):
    """Benjamini-Hochberg FDR q-values (monotone, clipped to [0,1])."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


# ---------- Fig. 1f: 8 parameters x 4 metrics, n = 9 ----------
df = pd.read_csv(os.path.join(DATA, 'Fig1f_plot_data_n9.csv'))
r_mat = pd.DataFrame(index=PARAMS, columns=LABELS, dtype=float)
p_mat = pd.DataFrame(index=PARAMS, columns=LABELS, dtype=float)
cells, all_p = [], []
for prm in PARAMS:
    for mi, met in enumerate(METRICS):
        r, p = stats.pearsonr(df[prm], df[met])
        r_mat.iloc[PARAMS.index(prm), mi] = r
        p_mat.iloc[PARAMS.index(prm), mi] = p
        cells.append((prm, LABELS[mi]))
        all_p.append(p)

q_all = bh(all_p)
q_mat = p_mat.copy()
for (prm, lab), q in zip(cells, q_all):
    q_mat.loc[prm, lab] = q

r_mat.to_csv(os.path.join(OUT, 'Fig1f_r_matrix.csv'))
p_mat.to_csv(os.path.join(OUT, 'Fig1f_p_matrix.csv'))
q_mat.to_csv(os.path.join(OUT, 'Fig1f_q_matrix.csv'))

print('=== Fig. 1f (n = %d) ===' % len(df))
print('tests:', len(all_p))
print('uncorrected P < 0.05:', int((p_mat.values < 0.05).sum()))
print('minimum BH q = %.4f' % q_all.min())
print()
print('Pearson r:')
print(r_mat.round(3).to_string())
print()
print('uncorrected P:')
print(p_mat.round(3).to_string())
print()
print('BH q:')
print(q_mat.round(3).to_string())

# ---------- Fig. 1c: onset concentration vs length, n = 13 ----------
d13 = pd.read_csv(os.path.join(DATA, 'transcripts_n13.csv'))
r_mass, p_mass = stats.pearsonr(d13['length'], d13['onset_conc_ng_uL'])
mol = d13['onset_conc_ng_uL'] / (330.0 * d13['length'])  # nM, 330 Da/nt
r_mol, p_mol = stats.pearsonr(d13['length'], mol)
print()
print('=== Fig. 1c (n = %d) ===' % len(d13))
print('mass onset conc. vs length : r = %+.2f, P = %.3f' % (r_mass, p_mass))
print('molar onset conc. vs length: r = %+.2f, P = %.3f' % (r_mol, p_mol))
print('onset conc. range: %.1f - %.1f ng/uL (fold %.2f)' % (
    d13['onset_conc_ng_uL'].min(), d13['onset_conc_ng_uL'].max(),
    d13['onset_conc_ng_uL'].max() / d13['onset_conc_ng_uL'].min()))
