# Internal analysis plan v5: revision analyses

Manuscript: bioengineering-4685597, "Robustness of a Machine-Learning-Derived Circulating MicroRNA Signature for Breast Cancer" (decision: reject, resubmission encouraged).
Written: 21 November 2026, before any analysis below was run. None of the statistics defined below had been computed when this file was hashed.
Status: internal plan. Its SHA-256 is recorded in `plan_v5_sha256.txt`.

## Inputs (frozen; SHA-256 recorded at the end of this file)

- Main analysis matrices (miRNA names on miRBase v22): `counts_matrix_umi_mirbase22.tsv`, `counts_matrix_reads_mirbase22.tsv`, and `qc_full.tsv` from the 19 October 2026 run. These produced the AUCs reported in the manuscript.
- Reproduction-run matrices (miRBase 22.1 from RNAcentral): `counts_matrix_umi.tsv`, `counts_matrix_reads.tsv`, `qc_full.tsv` from the 24 October 2026 run. These produced the nested cross-validation result (AUC 0.999).
- Locked signature: `locked_signature_v1.csv` (19 October 2026). Panels: Up8, Down11, LASSO5, High3 (let-7g-5p, miR-24-3p, miR-142-3p).
- Groups: 152 cases (BC001–BC152) and 48 female controls (HC001–HC048).

## R1. Permutation test of the nested cross-validated model (Reviewer 2, point 6)

Model: identical to plan v4, analysis B1. Features are miRNAs with at least 10 UMI molecules in at least 20% of libraries (a label-free filter), as log2(CPM + 1); standardised within training folds; L1-penalised logistic regression (liblinear, balanced class weights); C chosen from 12 log-spaced values from 0.01 to 10 by inner stratified 5-fold CV (random_state 100); performance from outer stratified 10-fold CV (random_state 0); statistic = pooled outer AUC.

1. Reproduction check: rerun B1 with 10 outer seeds. It must reproduce the recorded result (290 features; mean pooled AUC 0.999). If it does not, stop and report the discrepancy before step 2.
2. Observed statistic: pooled outer AUC with outer seed 0.
3. Null: 200 permutations of the case/control labels (numpy default_rng seed 2026, 152/48 preserved). Each permuted label vector goes through the full nested procedure of step 2, unchanged.
4. p = (1 + number of null AUCs ≥ observed) / 201. Report the null mean, SD, 95th percentile and maximum.
5. Decision rule: p < 0.05 means the near-perfect separation is not produced by the modelling pipeline itself (leakage or chance). Stated in advance: this test cannot say whether the separation is technical or biological.
6. Exploratory, declared in advance: B1 refitted among quality-matched libraries only (fraction of reads without adapter within the control range, 3.49–6.38%: 36 cases, 48 controls), same nested procedure, outer seed 0; reported without a decision rule.

## R2. Normalisation sensitivity of the panel scores (Reviewer 2, point 5)

Matrices: primary = main analysis matrices (UMI molecules and aligned reads); secondary = reproduction-run matrices.

Normalisations, each followed by log2(x + 1):
(a) CPM on total miRNA counts (as in the manuscript);
(b) DESeq2 median-of-ratios size factors (geometric means over miRNAs with non-zero counts in all 200 libraries);
(c) TMM (edgeR defaults: M-value trim 30%, A-value trim 5%, precision weights; reference library = the one whose upper-quartile-scaled counts are closest to the mean upper quartile);
(d) upper-quartile scaling (75th percentile of counts over miRNAs with at least one count in any library);
(e) quantile normalisation of log2(CPM + 1).

Scores, as in plan v4: (i) mean z-score of panel members (this study's score; z computed across all 200 libraries per miRNA; down members subtracted); (ii) mean log2 normalised expression (original authors' style). Panels: Up8, Down11 (members detected in the matrix), LASSO5, High3. AUC for 152 cases vs 48 controls with DeLong 95% CIs.

Decision rule for the down-regulated panel (same thresholds as plan v4): if the Down11 AUC reaches 0.75 or more in any normalisation × unit × score on the primary matrices, the panel's failure is "pipeline-sensitive"; if all are below 0.70, the failure is "robust to normalisation"; otherwise "partial". All grid values are reported, not only the best.

## R3. Alignment sensitivity (Reviewer 2, point 5; needs raw reads, run on a rented CPU server)

The 200 libraries are re-quantified with the same trimming and UMI extraction as the main pipeline and three alignment settings against the miRBase v22 mature sequences: (A) Bowtie 1, at most 1 mismatch (the manuscript's setting; acts as a check that the rerun reproduces the main matrix); (B) Bowtie 1, exact matches only; (C) Bowtie 2 end-to-end, `--very-sensitive`, `-L 10`, `--norc`. UMI molecules and aligned reads are counted for each. Each matrix is scored with normalisation (a) and both scores of R2. Same decision rule as R2, applied over settings B and C. If setting A does not reproduce the main-analysis Down11 AUC within ±0.02, this is reported and settings B and C are compared with A rather than with the manuscript.

## R4. Library quality by group (Reviewer 2, point 10)

Descriptive table from the main analysis `qc_full.tsv`: median and interquartile range by group of the eight library-quality metrics used by the technical-only classifier (raw reads; fraction without adapter; fraction with insert outside 16–30 nt; retained inserts; alignment rate; duplication rate; UMI molecules; miRNAs detected), with Mann–Whitney U p values, plus the classifier AUCs already reported in Section 3.4 (quality only, LASSO5 only, both, permuted-label null). No decision rule.

## R5. Cases outside the original discovery set (Reviewer 2, point 2)

Run only if the original authors supply the assignment of BC001–BC152 to their discovery and validation cohorts. Female-only AUCs (Up8, Down11, LASSO5, High3; score (i) on the main UMI matrix, normalisation (a)) for 76 validation cases vs 48 controls and for 76 discovery cases vs 48 controls, with DeLong 95% CIs; the difference (discovery − validation) with a 2,000-resample stratified bootstrap 95% CI. No decision rule; reported as an estimate of optimism. If no assignment is supplied, the manuscript states that it was requested and not available.

## What will not change

The locked signature, panel membership, the female-only design, and all results already reported. New results are added to the manuscript as revision analyses and labelled with this plan's date.

## Input SHA-256

```
d3f5680f91ea1d338cbc9de09ed6f8b05e714cffe7392faeaf4b213ea359e3f3  main/results/counts_matrix_umi_mirbase22.tsv
b41a2483448cdfd27cc47668fcb4acff4b3fc628886b37cd5d28b3dbfe9cbed2  main/results/counts_matrix_reads_mirbase22.tsv
64cf9407a930f9b5ab6b14e9e8e1a0a0b49c754728eb111ac783540a9ca87e07  main/results/qc_full.tsv
7970d078e3d0e5c92e5827b2ff1278a224e97e7123f16df8ad44e41bf8ae249d  main/results/locked_signature_v1.csv
e3204f2c62d17c430979649f5284415ba00a65a1c96a4ad144c9b51f69cfc3cd  rerun/counts_matrix_umi.tsv
5d8d72f145da7d316c73cf15dc537ccaef06dee295adae5b51de9b4b19a4d672  rerun/counts_matrix_reads.tsv
f3fab52fea32fb6fee053a8f5bacfcf2a28a602b299642839c0022a570e4a18e  rerun/results/qc_full.tsv
```
