#!/usr/bin/env python3
"""Run on the pod after r3_pod.sh matrix: keep only what scoring needs (panel miRNAs and library totals for all six
matrices; the 60 most abundant miRNAs of setting A UMI for the reproduction check; read-count completeness)."""
import pandas as pd, urllib.request, io
GH = "https://raw.githubusercontent.com/MohannadBanat/mena-mirna-signature-robustness/main"
sig = pd.read_csv(io.StringIO(urllib.request.urlopen(f"{GH}/analysis_plan/locked_signature_v1.csv").read().decode()))
panel = sorted(set(sig["miRNA ID"]) | {"hsa-let-7g-5p", "hsa-miR-24-3p", "hsa-miR-142-3p"})
out = []
for X in "ABC":
    for u in ["umi", "reads"]:
        m = pd.read_csv(f"r3_{X}_{u}.tsv", sep="\t", index_col=0)
        keep = [p for p in panel if p in m.index]
        sub = m.loc[keep].copy(); sub.loc["__TOTAL__"] = m.sum(); sub.loc["__DETECTED__"] = (m > 0).sum()
        if X == "A" and u == "umi":
            top = m.div(m.sum()).mean(1).sort_values(ascending=False).index[:60]
            m.loc[top].to_csv("sub_A_umi_top60.tsv", sep="\t")
        sub.to_csv(f"sub_{X}_{u}.tsv", sep="\t")
runs = pd.read_csv("runs.tsv", sep="\t"); qc = pd.read_csv("r3_qc.tsv", sep="\t")
mm = qc.merge(runs, left_on="run", right_on="run_accession")
mm[["sample", "run", "reads", "read_count", "kept", "no_adapter", "bad_len"]].to_csv("sub_qc.tsv", sep="\t", index=False)
print("libraries:", len(mm), "read counts equal to ENA:", int((mm.reads == mm.read_count).sum()))
