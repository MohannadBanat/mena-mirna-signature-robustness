# Amendment v4a to plan v4 — written 2026-10-24, before the rerun produced any data
1. Correction: the UMI panel scores in analysis.py use log2(CPM + 1) (counts per million molecules), not size-factor-normalised
   counts as plan v4 and Methods 2.5 state. A1/A2 therefore use log2(CPM + 1) for both read and UMI counts, exactly as analysis.py.
2. Added A5 (secondary): original-style panel scores (mean log2(CPM + 1) without z-scoring) from read counts and from UMI counts,
   to separate the effect of the counting unit from the effect of the scoring method.
3. Cost cap: the pod stops itself 7 hours after the job starts.
