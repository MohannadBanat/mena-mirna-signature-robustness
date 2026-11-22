# Amendment v5a to internal analysis plan v5

Written 22 November 2026, after the plan v5 results R1, R2 and R4 were known (results_v5_summary.md, 21 November 2026) and before any analysis below was run. Plan v5 itself is unchanged; its SHA-256 is in plan_v5_sha256.txt.

## Why this amendment

R2 showed that the panel AUCs depend on normalisation (Down11 0.66 with counts per million, 0.78 with DESeq2 median-of-ratios). Three further results in the manuscript used counts-per-million scores. The authors decided on 22 November 2026 to report counts per million (the normalisation fixed in the 19 October 2026 plan) and DESeq2 median-of-ratios (the original study's normalisation) side by side, without changing the primary analysis.

## A1. Side-by-side recomputation (descriptive; no decision rules)

Score: LASSO5 = mean z-score (across all 200 libraries) of the five signature miRNAs, computed from log2(normalised count + 1) on the main-analysis UMI matrix. Normalisations: (a) counts per million (as reported; must reproduce the published per-sample scores within rounding, else stop and report) and (b) DESeq2 median-of-ratios exactly as in plan v5, R2 (b).

1. LASSO5 among quality-matched libraries (fraction of reads without adapter within the control range, 3.49–6.38%; 36 cases, 48 controls): AUC with a 2,000-resample bootstrap 95% CI (numpy default_rng seed 0), as in qc_negative_control.py.
2. Technical-only comparison, as in qc_negative_control.py: logistic regression on the LASSO5 score alone and on the eight library-quality metrics plus the LASSO5 score; 20 repeats of stratified 10-fold CV (random_state 0–19); mean and range of the CV AUC.
3. BC144–BC152 run, as described in Section 3.4: mean score of BC144–BC152, of the other 143 cases and of the 48 controls; percentage of cases and of controls scoring below zero; the probability that nine randomly drawn cases have a mean score at or below the BC144–BC152 mean (100,000 draws without replacement, seed 0); AUC excluding BC144–BC152 with a 2,000-resample bootstrap 95% CI (seed 0); Spearman correlation between score and UMI molecules among cases.
4. Figure 3 and Table 3 redrawn with both normalisations.

Counts per million remains the primary analysis. DESeq2 median-of-ratios is labelled as the like-for-like comparison with the original study.

## A2. Implementation of plan v5, R3 (alignment sensitivity)

Run on a rented CPU server. Trimming, UMI capture and insert-length filter (16–30 nt) exactly as in pipeline/pod_pipeline.sh. Reference: human mature miRNA sequences from RNAcentral's miRBase mapping, built with rerun_2026-10-24/mkref.py; the number of sequences and the SHA-256 of the FASTA are recorded. Settings, all with --norc and one reported alignment per read:
(A) bowtie 1.3.1 -n 1 -l 16 --best --strata -k 1 (the pipeline setting);
(B) bowtie 1.3.1 -v 0 --best --strata -k 1 (exact matches only);
(C) bowtie2 end-to-end --very-sensitive -L 10 -N 0 (best alignment; ties broken by bowtie2's default random choice, seed 0).
UMI molecules by umi_tools 1.1.5 dedup --per-contig --method=unique; aligned reads from samtools idxstats. In addition to the plan v5 check against the main analysis, setting A is compared with the 24 October 2026 reproduction-run matrix (per-library Spearman correlation of counts and the Down11 AUC). Scoring and the decision rule are as in plan v5, R3.
