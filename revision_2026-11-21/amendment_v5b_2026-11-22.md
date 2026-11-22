# Amendment v5b to internal analysis plan v5 (erratum to amendment v5a, A2)

Written 22 November 2026, during setup of the alignment analysis (plan v5, R3) and before any alignment result was produced.

Amendment v5a, A2, described setting (A) as "bowtie 1.3.1 -n 1 -l 16 --best --strata -k 1 (the pipeline setting)" and the UMI step as "umi_tools dedup --per-contig". Those flags were copied from an outdated draft of the pipeline. The pipeline that produced the manuscript's matrices (pipeline/pod_pipeline.sh and rerun_2026-10-24/pod_pipeline.sh in this repository) uses:

- alignment: bowtie -v 1 --best -k 1 --norc (end-to-end, at most one mismatch);
- UMI collapsing: umi_tools dedup --per-gene --per-contig --method=unique --umi-separator=_.

Setting (A) is therefore run with the committed pipeline's flags, as A2 intended ("the pipeline setting"). Setting (B) becomes bowtie -v 0 --best -k 1 --norc (A's flags with zero mismatches). Setting (C) and everything else in A2 are unchanged. --strata is not used (it has no effect with -k 1).

Data source: ENA served the server at about 0.15 MB/s on 22 November 2026, so reads were taken from NCBI's SRA mirror on AWS (sra-pub-run-odp) and converted with fasterq-dump 3.1.1 (--split-spot). Completeness is checked against ENA's read counts; the ENA MD5 checksums do not apply to these files. For the test library (BC090), the read count (11,417,414) and the trimming counts (3,813,820 retained; 1,128,520 without adapter; 6,475,074 outside 16-30 nt) were identical to the 24 October 2026 run.
