# Internal analysis plan v4 — written 2026-10-24, before the rerun produced any data
Pod: RunPod "mRNA rerun" (vtz6ktw4i3sz7s), deployed 2026-10-24. Reprocessing of PRJNA1436574 with the unchanged pipeline
(pod_pipeline.sh: QIAseq adapter, 12-nt UMI, 16-30 nt inserts, Bowtie 1.3.1 -v 1 --best -k 1 --norc, UMI-tools 1.1.5 dedup per miRNA,
MIMAT mapping to miRBase v22 names). Only change: samples are processed in parallel. The rerun must first reproduce the published
re-analysis (1,864 miRNAs; 37/38 replicated; LASSO5 AUC 0.92; down panel AUC 0.66) before the new analyses are interpreted.

## A. Read counts versus UMI molecules (tests the Section 4.2 hypothesis)
Read counts = aligned reads per miRNA before UMI collapsing (same libraries, same alignment).
A1 (primary). Down-regulated panel score (negative mean z-score of the 9 detectable panel miRNAs, log2 size-factor-normalised count + 1,
    exactly as for UMI counts) for 152 cases vs 48 controls, computed from read counts.
    Decision: hypothesis SUPPORTED if read-based AUC >= 0.75 (closer to the published 0.85 than to our UMI 0.66);
    NOT SUPPORTED if read-based AUC < 0.70; otherwise PARTIAL.
A2. Same for the up-regulated panel, LASSO5 and the low-count LASSO pair (miR-138-5p, miR-32-3p).
A3. Female-only DESeq2 on read counts with the same filter: number of the 19 panel miRNAs that pass the filter and replicate,
    vs UMI counts; Spearman rho of log2FC with published log2FC, all shared miRNAs and published baseMean < 10.
A4. Reads per molecule (duplication) per miRNA: median across libraries by abundance class, and case/control difference
    (Mann-Whitney) for the low-count panel miRNAs.

## B. Nested cross-validated LASSO re-derivation (what a properly validated model achieves)
Features: log2(CPM + 1) per library (no cross-sample information) of miRNAs with >= 10 molecules in >= 20% of libraries (label-free filter).
Model: L1-penalised logistic regression, features standardised within training folds, class weights balanced; C chosen by inner
stratified 5-fold CV (AUC) over 12 log-spaced values 0.01-10; outer stratified 10-fold CV repeated with 10 seeds. Report mean outer AUC
and range, and selection frequency of the five LASSO5 miRNAs. Variants: B1 UMI counts (primary); B2 read counts; B3 UMI miRNA features
plus the eight library-quality metrics. Benchmark: the technical-only classifier (cross-validated AUC 0.94).
B is descriptive; no pass/fail criterion.

## C. GSE73002 re-check
Rerun gse.py unchanged to confirm Table 4 and the PPV range (0.5% prevalence, control specificity 3.8-4.4%).
