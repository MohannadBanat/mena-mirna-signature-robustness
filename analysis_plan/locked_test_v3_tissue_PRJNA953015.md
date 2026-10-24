# Locked plan v3: MENA tumour-tissue check of the plasma signature miRNAs (PRJNA953015)
Written before any read from PRJNA953015 was downloaded. Seen so far: ENA run metadata only
(96 single-end miRNA-Seq runs, HiSeq 4000, sample titles miRNA_BC##; reads per run 1.6M-20.7M).
Cohort (Elango et al., Cell Death Dis 2023;14:415): 96 FFPE treatment-naive invasive breast tumours,
Hamad Medical Corporation, Qatar; 60 MENA and 36 non-MENA patients; QIAseq miRNA library kit.
Note: the paper lists PRJNA954402 for miRNA data, but the deposit shows PRJNA953015 is the miRNA-Seq project.

## Processing (identical to the plasma pipeline)
Adapter AACTGTAGGCACCATCAAT, 12-nt UMI, inserts 16-30 nt, Bowtie 1.3.1 (-v 1 --best -k 1 --norc) to the same miRBase
reference, UMI-tools dedup per miRNA, MIMAT mapping to miRBase v22 names. Molecules per million (MPM) = UMI molecules / total
molecules x 1e6.

## Questions and pre-specified readouts (descriptive; no tissue-vs-normal comparison is possible)
Q1 Presence in tumours. For each of the 19 panel miRNAs: % of tumours with >=1 molecule, % with >=10 MPM, median MPM and
   abundance rank among all miRNAs. A plasma-up miRNA is called "supported as tumour-expressed" if >=10 MPM in >=50% of tumours,
   "rare in tumours" if detected (>=1 molecule) in <50% of tumours, otherwise "low".
Q2 miR-4455 specifically: classified by the Q1 rule; this is the decisive readout for the plasma all-or-none pattern.
Q3 Immune origin. Spearman correlation across tumours of log2(MPM+1) of miR-142-3p and miR-223-3p with miR-150-5p
   (lymphocyte-enriched reference); rho >= 0.5 is read as consistent with an immune-cell source.
Q4 Plasma direction vs tissue abundance. Among miRNAs tested in the female-only plasma analysis, compare tumour median MPM of
   plasma-up versus plasma-down DE miRNAs (Mann-Whitney); exploratory.
Q5 MENA restriction. If per-sample ethnicity is available, Q1-Q2 are repeated in the 60 MENA tumours; if not, all 96 are used and
   this is stated.
No threshold, miRNA list or rule will be changed after the data are seen; any further analysis is labelled exploratory.
