# Revision analyses (21–22 November 2026)

Analyses added after peer review of manuscript bioengineering-4685597. They use the matrices and locked signature of the main analysis and change no result reported in the earlier version.

## Plans (written and pushed before the analyses they describe were run)
| File | Written | Commit | SHA-256 file |
|---|---|---|---|
| `analysis_plan_v5_revision.md` | 21 Nov 2026 | 0c666d5 (pushed 18:29 UTC) | `plan_v5_sha256.txt` |
| `amendment_v5a_2026-11-22.md` (side-by-side normalisations; alignment run) | 22 Nov 2026 | 6fec93e | `amendment_v5a_sha256.txt` |
| `amendment_v5b_2026-11-22.md` (erratum: alignment flags of the committed pipeline) | 22 Nov 2026 | 8c7c222 | `amendment_v5b_sha256.txt` |

Note on v5b: it says reads were converted with `fasterq-dump --split-spot`. During the run, fasterq-dump's temporary files filled the disk, so reads were streamed with `fastq-dump --stdout` instead (commit e12f13c, `r3_pod.sh`). Completeness was checked against ENA's read counts for all 200 libraries.

## Scripts and outputs
| Analysis (plan item) | Script | Output in `out/` |
|---|---|---|
| R1 label-permutation test of the nested model | `r1_permutation.py` | `r1_*.json`, `r1_null.txt` |
| R2 normalisation sensitivity (CPM, median-of-ratios, TMM, upper-quartile, quantile) | `r2_normalisation.py` | `r2_normalisation_grid.tsv`, `r2_decision.json` |
| R3 alignment sensitivity (Bowtie 1 ≤1 mismatch, Bowtie 1 exact, Bowtie 2 end-to-end) | `r3_pod.sh` (run on a RunPod CPU server), `r3_pod_extract.py`, `r3_score.py` | `r3_alignment_grid.tsv`, `r3_decision.json`; matrices in `r3_matrices/` |
| R4 single library-quality metrics | `r4_quality_table.py` | `r4_metric_aucs.tsv` |
| A1 LASSO5 with CPM and median-of-ratios side by side | `a1_side_by_side.py` | `a1_side_by_side.json`, `a1_lasso5_scores_*.tsv` |
| Figure 3 (revised) | `fig3_revised.py` | `Figure3_revised.png` |
| Summary | `summarise.py` | `results_v5_summary.md` |

R5 (AUCs for the 76 cases outside the original discovery set) requires the original authors' cohort assignment, which was requested and was not available at resubmission.

`r3_matrices/` holds the six count matrices of the alignment analysis (settings A–C × UMI molecules and aligned reads; 200 libraries each), the per-library trimming counts (`r3_qc.tsv`), the reference size and SHA-256 (`reference_info.txt`) and tool versions. `r3_results_sha256_on_pod.txt` is the checksum of the results archive computed on the server.

Scripts use absolute paths from the analysis machine (`/home/claude/rev/...`); set `W` at the top of each script to your working directory. Inputs are the main-analysis matrices regenerated with `../pipeline/` and the reproduction-run matrices from `../rerun_2026-10-24/`; their SHA-256 checksums are listed in the plan.
