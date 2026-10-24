# Locked test v2: benign breast disease specificity of the MENA signature in two further cohorts
Written before any expression data or sample labels from GSE255660 or GSE225117 were downloaded.
Only the GEO series titles, summaries and sample-title counts (207 "Breast Benign Disease", 58 "Breast Cancer",
30 "non-cancer healthy donor", 24 technical replicates for GSE255660) had been seen.

## Signature (unchanged from locked_signature_v1.csv)
LASSO5 (all up in cancer): hsa-let-7g-5p, hsa-miR-138-5p, hsa-miR-142-3p, hsa-miR-24-3p, hsa-miR-32-3p
High3: hsa-let-7g-5p, hsa-miR-24-3p, hsa-miR-142-3p
Up8: LASSO5 + hsa-miR-223-3p, hsa-miR-4455, hsa-miR-766-5p

## GSE255660 (serum, 3D-Gene V22, women with abnormal mammograms; Akiyama et al. 2025)
Samples: primary analysis uses one array per subject: "Breast Cancer" (cases) vs "Breast Benign Disease" (negatives).
Healthy-donor and technical-replicate arrays are excluded from the primary analysis.
If annotation separates DCIS from invasive cancer: primary = all cancers as labelled; secondary = invasive only.
Expression: series-matrix normalised values, log2-transformed if not already on log scale; probes matched by MIMAT accession or miRBase name.
Scores:
 S1 (primary) LASSO5 direction-only: mean over available LASSO5 miRNAs of z = (x - mean_benign) / sd_benign. No fitted weights.
 S2 High3 direction-only (same construction).
 S3 Up8 direction-only (miRNAs available on the array).
 S4 Best case for the signature: logistic regression on the LASSO5 miRNAs trained on cancer vs benign within GSE255660,
    5-fold stratified cross-validation repeated 20 times; out-of-fold predictions averaged across repeats.
Positive control (does the cohort contain any cancer-vs-benign signal?):
 PC L1-penalised logistic regression on all miRNAs, nested 5-fold CV (inner 5-fold CV over C), 10 repeats.
Failure criterion (as for GSE73002): the signature fails if, for S1, the 95% DeLong CI of the AUC includes 0.5
 OR specificity at 90% sensitivity is below 50%. S2-S4 are reported regardless.
PPV at 90% sensitivity at 10% (symptomatic) and 0.5% (screening) prevalence.
Descriptive only: S1 in the 5 healthy donors; S1 reproducibility across technical replicates.

## GSE225117 (plasma, small RNA sequencing; Khadka et al. 2023)
Groups: early-stage cancer (CA) vs all benign (high-, moderate- and no-risk benign) pooled.
Expression: processed counts as deposited, normalised to counts per million, log2(CPM + 1).
Scores S1-S3 as above (z relative to the pooled benign group). Small n: AUC with DeLong CI reported as supportive only;
the failure criterion is applied but a pass or fail here does not override GSE255660.

## Not changed
No miRNA, weight, threshold or criterion will be changed after the data are seen. Any additional analysis will be labelled exploratory.
