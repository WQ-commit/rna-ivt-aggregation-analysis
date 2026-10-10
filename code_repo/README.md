# Correlation and Benjamini–Hochberg analysis

Analysis code for the correlation analyses and multiple-testing correction
reported in **Fig. 1c** and **Fig. 1f** of the manuscript.

## Contents

```
correlation_bh_analysis.py   main script (Pearson r, two-sided P, BH q)
data/Fig1f_plot_data_n9.csv  9 diverse transcripts, 8 parameters x 4 metrics
data/transcripts_n13.csv     13 transcripts, onset concentration vs length
output/                      generated result matrices (r, P, q)
```

## What it computes

- **Fig. 1f** — Pearson correlation between 8 sequence/structure parameters
  (length, GC, MFE/nt, paired fraction, stems, max loop, LZ76, k2_H) and
  4 aggregation metrics (apparent onset time, apparent onset concentration,
  mass-unit rate, molar rate) across the 9 diverse transcripts, i.e.
  8 x 4 = 32 tests. Benjamini–Hochberg FDR correction is applied across all
  32 tests; three cells reach uncorrected P < 0.05 and none survives
  correction.
- **Fig. 1c** — mass onset concentration vs length and molar onset
  concentration vs length across the 13 transcripts.

## Requirements

Python 3.8+ with `numpy`, `pandas`, `scipy`.

```
pip install -r requirements.txt
```

## Run

```
python correlation_bh_analysis.py
```

Results are written to `output/` and printed to the console.

## Data

- `Fig1f_plot_data_n9.csv` — the 9-transcript set used for the correlation
  matrix (Fig. 1f). Onset concentrations and rates are technical-replicate
  means (n = 3).
- `transcripts_n13.csv` — the full 13-transcript set with apparent onset
  time, apparent onset concentration and 150-min final concentration.

## Citation

If you use this code, please cite the manuscript (reference to be added
upon publication).

## Repository

- GitHub: `<ADD GITHUB URL HERE>`
- Zenodo DOI: `<ADD ZENODO DOI HERE>`
