# Robustness of a machine-learning-derived circulating microRNA signature for breast cancer — code and analysis plans

Code accompanying: Abu Owida H., Alsarhan A., Banat M., Al-Fraihat D., Algarni A. "Robustness of a Machine-Learning-Derived Circulating MicroRNA Signature for Breast Cancer: Confounding, Library Quality and Benign-Disease Specificity." Submitted to *Bioengineering* (MDPI), 2026.

## Contents
- `pipeline/` — download of PRJNA1436574 (ENA), read processing (QIAseq adapter trimming, 12-nt UMI extraction, 16–30 nt inserts, Bowtie 1.3.1, UMI-tools 1.1.5 dedup, miRBase mapping), and the same pipeline for the tumour cohort PRJNA953015 (`t.sh`).
- `analysis/` — female-only differential expression (PyDESeq2), depth matching by binomial thinning, haemolysis reassessment and panel scores (`analysis.py`); GSE73002 benign-disease test (`gse.py`, `roc.py`); tumour-tissue readouts (`ta.py`, `tissue_an.py`); QC table (`qctab.py`); exploratory technical-only negative-control classifier (`qc_negative_control.py`, added 24 October 2026 after the main results were known; reads `data/qc_full.txt` and `data/scores.txt`).
- `figures/` — scripts for Figures 1–6.
- `analysis_plan/` — internal analysis plans. `bundle.tgz` (19 October 2026, 00:39 UTC) holds the analysis scripts and `locked_signature_v1.csv` as committed before they were run; `locked_test_v3_tissue_PRJNA953015.md` (21 October 2026, 16:36:58 UTC) is the tumour-tissue plan, written before any tumour read was downloaded. SHA-256 hashes are in `lock_timestamp.txt` and `bundle_sha256.txt`. `locked_test_v2_GSE255660_GSE225117.md` and `amendment_1.md` record a planned test on two non-MENA cohorts that was not run after the study was restricted to MENA cohorts plus GSE73002. These plans were not publicly registered.
- `data/` — derived summary tables used in the manuscript. `pub_s1.tsv` and `pub_s15.tsv` are transcribed from Supplementary Tables S1 and S15 of Vishnubalaji et al., J. Transl. Med. 2026, 24, 1233 (doi:10.1186/s12967-026-08452-z), reproduced under that article's CC BY licence.

## Raw data
- Plasma: BioProject PRJNA1436574 (SRA/ENA).
- Benign breast disease: GEO GSE73002.
- Tumours: BioProject PRJNA953015.

The UMI count matrix (1,864 miRNAs × 200 libraries) is regenerated from the raw reads with `pipeline/`. Software: Python 3.11, NumPy, SciPy, pandas, scikit-learn, PyDESeq2 0.5.4.

## Reproduction run (24 October 2026)
`rerun_2026-10-24/` holds plan v4 and amendment v4a (with SHA-256 hashes, written before the rerun produced data), the orchestrator
(`run_all.sh`), the pipeline as run (with a fallback that builds the miRBase 22.1 human reference from RNAcentral, used because
mirbase.org was unreachable that day; `mkref.py`), the read-count versus UMI and nested cross-validation analyses (`extra.py`), and
the key outputs (`rerun_results_2026-10-24.txt`). Under plan v4 the UMI-collapsing explanation for the down-regulated panel was not supported.
