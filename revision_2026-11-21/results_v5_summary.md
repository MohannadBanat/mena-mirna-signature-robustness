# Revision analyses (plan v5, pushed to GitHub 2026-11-21 18:29 UTC before these were run)

## R1 Nested model permutation test
- Reproduction check (10 seeds): 290 features, mean AUC 0.9989 (range 0.9984-0.9993); recorded 0.999 (0.998-0.999): reproduced.
- Observed AUC (outer seed 0): 0.9992
- Null over 200 label permutations: mean 0.4965, SD 0.0649, 95th percentile 0.603, max 0.6796
- Null AUCs >= observed: 0; p = 0.00498 (decision rule p < 0.05: met)
- Exploratory, quality-matched libraries (no-adapter fraction 3.49-6.38%): 36 cases, 48 controls, nested AUC 1.0

## R2 Normalisation sensitivity
- Down11 on main matrices (UMI and reads, both scores, CPM/median-of-ratios/TMM/quantile): 0.644 to 0.837 -> pre-specified verdict: pipeline-sensitive (>= 0.75 reached with median-of-ratios, TMM and quantile)
- Upper-quartile scaling: undefined in 36 of 200 libraries (75th percentile of counts = 0); its cells are degenerate and excluded from interpretation.

Main matrices, UMI molecules, z-score (AUC, 95% CI):

| Normalisation | Up8 | Down11 | LASSO5 | High3 |
|---|---|---|---|---|
| CPM | 0.925 (0.887-0.963) | 0.660 (0.575-0.746) | 0.922 (0.882-0.962) | 0.890 (0.836-0.944) |
| median-of-ratios | 0.958 (0.931-0.985) | 0.779 (0.703-0.854) | 0.957 (0.927-0.986) | 0.935 (0.896-0.975) |
| TMM | 0.963 (0.937-0.989) | 0.740 (0.660-0.819) | 0.962 (0.936-0.988) | 0.936 (0.898-0.974) |
| quantile | 0.964 (0.939-0.988) | 0.807 (0.735-0.879) | 0.964 (0.941-0.987) | 0.947 (0.913-0.981) |

Exploratory diagnostic (added after R2, not in plan v5): miR-16-5p makes up about half of the molecules in a typical library (median 51.5%, Table 1), so CPM values of every other miRNA move with its share. Median-of-ratios size factors relative to library size differed between groups (median log2 ratio 0.39 in cases vs 0.08 in controls, Mann-Whitney p = 0.0003), i.e. CPM and composition-robust normalisations disagree systematically by group.

## R4 Single library-quality metrics, cases vs controls (AUC in either direction, DeLong 95% CI)

| Metric | Direction | AUC (95% CI) |
|---|---|---|
| Raw reads | higher in cases | 0.735 (0.651-0.819) |
| Fraction without adapter | higher in cases | 0.925 (0.890-0.960) |
| Fraction with insert outside 16-30 nt | higher in cases | 0.526 (0.420-0.633) |
| Retained inserts | higher in cases | 0.641 (0.548-0.734) |
| Alignment rate | lower in cases | 0.645 (0.560-0.730) |
| Duplication rate | lower in cases | 0.501 (0.410-0.591) |
| UMI molecules | higher in cases | 0.595 (0.504-0.686) |
| miRNAs detected | higher in cases | 0.655 (0.567-0.742) |

## A1 Side by side: counts per million vs DESeq2 median-of-ratios (amendment v5a, pushed 2026-11-22 before running)

CPM column reproduces the published values exactly (per-sample scores within 0.0005).

| Result | Counts per million (published) | DESeq2 median-of-ratios |
|---|---|---|
| LASSO5 AUC, all 200 (bootstrap 95% CI) | 0.922 (0.878-0.958) | 0.957 (0.924-0.983) |
| LASSO5 AUC, quality-matched (36 cases, 48 controls) | 0.932 (0.878-0.975) | 0.968 (0.932-0.995) |
| CV AUC, LASSO5 alone (mean, range over 20 x 10-fold) | 0.917 (0.914-0.920) | 0.954 (0.951-0.956) |
| CV AUC, library quality + LASSO5 | 0.975 (0.965-0.981) | 0.977 (0.969-0.981) |
| BC144-BC152 mean score (all nine below zero: True / True) | -0.74 | -0.45 |
| Other 143 cases, mean score | 0.33 | 0.31 |
| Controls, mean score | -0.84 | -0.83 |
| Random nine-case draws with mean <= cluster mean (of 100,000) | 0 | 0 |
| LASSO5 AUC excluding BC144-BC152 | 0.943 (0.904-0.974) | 0.968 (0.939-0.990) |
| Spearman score vs UMI molecules, cases (rho, p) | -0.26, 0.0011 | -0.27, 0.00096 |
| Cases / controls scoring below zero | 32% / 96% | 32% / 96% |

Error found in the submitted manuscript (Section 3.4): it says 'only 4% of controls scored below zero'; 46 of 48 controls (96%) scored below zero, so it should read 'only 4% of controls scored above zero'.

## R3 Alignment sensitivity (plan v5 R3; amendments v5a A2 and v5b; run 22 Nov 2026 on RunPod, 16 vCPU, $0.65/h)

- Reads: all 200 libraries, read counts identical to ENA's (taken from NCBI's SRA mirror on AWS; see amendment v5b). Reference: 2,656 human mature miRNAs (RNAcentral miRBase mapping), SHA-256 f4197e44...
- Setting A (pipeline flags) reproduces the 24 Oct reproduction run: Down11 AUC 0.633 vs 0.633; per-library Spearman with the 24 Oct UMI matrix median 1.0, min 0.9876 (1791 miRNAs). It differs from the manuscript's main-analysis value (0.660) by more than 0.02, as the 24 Oct run already did (different reference source), so, as pre-specified, B and C are compared with A.
- Down11 AUC across settings B (exact match) and C (Bowtie 2), both units, both scores, counts per million: 0.623 to 0.656 -> pre-specified verdict: robust to alignment.

| Setting | Unit | Up8 | Down11 | LASSO5 | High3 |
|---|---|---|---|---|---|
| A: Bowtie 1, <=1 mismatch (pipeline) | UMI molecules | 0.924 | 0.633 | 0.921 | 0.888 |
| A: Bowtie 1, <=1 mismatch (pipeline) | aligned reads | 0.919 | 0.632 | 0.914 | 0.891 |
| B: Bowtie 1, exact match | UMI molecules | 0.903 | 0.637 | 0.921 | 0.887 |
| B: Bowtie 1, exact match | aligned reads | 0.897 | 0.628 | 0.913 | 0.888 |
| C: Bowtie 2 end-to-end | UMI molecules | 0.925 | 0.626 | 0.922 | 0.891 |
| C: Bowtie 2 end-to-end | aligned reads | 0.919 | 0.623 | 0.916 | 0.893 |

(z-score, counts per million. Full grid with mean-log scores and DeLong CIs: out/r3_alignment_grid.tsv.)

Interpretation: the down-regulated panel's shortfall does not depend on alignment (0.62-0.66 across three aligners/settings) but does depend on normalisation (R2: 0.66 with CPM, 0.78 with DESeq2 median-of-ratios).
