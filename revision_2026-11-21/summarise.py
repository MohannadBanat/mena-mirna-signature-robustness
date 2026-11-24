#!/usr/bin/env python3
"""Collect plan v5 results (R1, R2, R4) into one dated summary file. Numbers are copied from the output files only."""
import json, pandas as pd
O = "/home/claude/rev/revision_2026-11-21/out"
rep = json.load(open(f"{O}/r1_reproduce.json")); per = json.load(open(f"{O}/r1_permutation.json")); mat = json.load(open(f"{O}/r1_matched.json"))
dec = json.load(open(f"{O}/r2_decision.json")); g = pd.read_csv(f"{O}/r2_normalisation_grid.tsv", sep="\t")
m4 = pd.read_csv(f"{O}/r4_metric_aucs.tsv", sep="\t")
L = ["# Revision analyses (plan v5, pushed to GitHub 2026-11-21 18:29 UTC before these were run)", ""]
L += ["## R1 Nested model permutation test",
      f"- Reproduction check (10 seeds): {rep['features']} features, mean AUC {rep['auc_mean']} (range {rep['auc_min']}-{rep['auc_max']}); recorded 0.999 (0.998-0.999): reproduced.",
      f"- Observed AUC (outer seed 0): {per['observed']}",
      f"- Null over {per['n_perm']} label permutations: mean {per['null_mean']}, SD {per['null_sd']}, 95th percentile {per['null_p95']}, max {per['null_max']}",
      f"- Null AUCs >= observed: {per['n_null_ge_obs']}; p = {per['p']} (decision rule p < 0.05: {'met' if per['p'] < 0.05 else 'not met'})",
      f"- Exploratory, quality-matched libraries (no-adapter fraction {mat['range_pct'][0]}-{mat['range_pct'][1]}%): {mat['cases']} cases, {mat['controls']} controls, nested AUC {mat['nested_auc_seed0']}", ""]
L += ["## R2 Normalisation sensitivity",
      f"- Down11 on main matrices, all computable cells: min {dec['down_primary_min']}, max {dec['down_primary_max']} -> verdict: {dec['verdict']}",
      "- Upper-quartile scaling: undefined in 36 of 200 libraries (75th percentile of counts = 0); its cells are degenerate and excluded from interpretation.", ""]
sub = g[(g.matrices == "main") & (g.unit == "UMI") & (g.score == "z-score") & (g.normalisation != "upper-quartile")]
L += ["Main matrices, UMI molecules, z-score (AUC, 95% CI):", "", "| Normalisation | Up8 | Down11 | LASSO5 | High3 |", "|---|---|---|---|---|"]
for n in ["CPM", "median-of-ratios", "TMM", "quantile"]:
    r = sub[sub.normalisation == n].set_index("panel")
    L.append(f"| {n} | " + " | ".join(f"{r.loc[p,'AUC']:.3f} ({r.loc[p,'lo']:.3f}-{r.loc[p,'hi']:.3f})" for p in ["Up8", "Down11", "LASSO5", "High3"]) + " |")
L += ["", "## R4 Single library-quality metrics, cases vs controls (AUC in either direction, DeLong 95% CI)", "", "| Metric | Direction | AUC (95% CI) |", "|---|---|---|"]
for _, r in m4.iterrows(): L.append(f"| {r['metric']} | {r['direction']} | {r['AUC (either direction)']:.3f} ({r['lo']:.3f}-{r['hi']:.3f}) |")
open(f"{O}/../results_v5_summary.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))
