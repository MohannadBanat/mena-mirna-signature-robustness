import pandas as pd, numpy as np
m=pd.read_csv("tissue_counts_umi_v22.tsv",sep="\t",index_col=0)
q=pd.read_csv("tissue_qc.tsv",sep="\t")
mpm=m.div(m.sum(),axis=1)*1e6; rk=mpm.median(1).rank(ascending=False,method="min")
P=["hsa-let-7g-5p","hsa-miR-138-5p","hsa-miR-142-3p","hsa-miR-24-3p","hsa-miR-32-3p","hsa-miR-223-3p","hsa-miR-4455","hsa-miR-766-5p","hsa-miR-106b-3p","hsa-miR-4497","hsa-miR-660-3p","hsa-miR-30a-5p","hsa-miR-3195","hsa-miR-500b-5p","hsa-miR-1976","hsa-miR-106a-3p","hsa-miR-382-3p","hsa-miR-4639-5p","hsa-miR-6741-3p"]
REF=["hsa-miR-21-5p","hsa-miR-150-5p","hsa-miR-155-5p","hsa-miR-451a","hsa-miR-23a-3p","hsa-miR-16-5p","hsa-miR-122-5p","hsa-miR-486-5p","hsa-miR-205-5p","hsa-miR-200c-3p"]
UP="hsa-miR-122-5p,hsa-miR-130a-3p,hsa-miR-130b-3p,hsa-miR-142-3p,hsa-miR-151a-5p,hsa-miR-151b,hsa-miR-22-5p,hsa-miR-221-3p,hsa-miR-24-3p,hsa-miR-28-5p,hsa-miR-320a-3p,hsa-miR-337-5p,hsa-miR-34a-5p,hsa-miR-369-3p,hsa-miR-4433b-5p,hsa-miR-4455,hsa-miR-4521,hsa-miR-487b-3p,hsa-miR-491-5p,hsa-miR-493-5p,hsa-miR-505-5p,hsa-miR-574-5p,hsa-miR-627-5p,hsa-miR-6734-5p,hsa-miR-874-3p".split(",")
DN="hsa-miR-1270,hsa-miR-146b-3p,hsa-miR-146b-5p,hsa-miR-181a-2-3p,hsa-miR-29a-3p,hsa-miR-30d-5p,hsa-miR-30e-3p,hsa-miR-3168,hsa-miR-342-3p,hsa-miR-532-5p,hsa-miR-660-5p,hsa-miR-9-5p".split(",")
o=open("ta_results.txt","w")
w=lambda *a: o.write("\t".join(str(x) for x in a)+"\n")
w("QC", m.shape[0], m.shape[1], int(q.reads.median()), round((q.kept/q.reads).median()*100,1), int(q.umi_total.median()), int(q.mirnas_detected.median()), int(q.umi_total.min()), int(q.umi_total.max()), int(q.reads.sum()), int(q.umi_total.sum()))
w("TOP15", ";".join(f"{i.replace('hsa-','')}:{mpm.loc[i].median():.0f}" for i in mpm.median(1).sort_values(ascending=False).index[:15]))
for g,L in [("PANEL",P),("REF",REF),("UP",UP),("DN",DN)]:
    for p in L:
        if p in m.index: w(g,p.replace("hsa-",""),f"{(m.loc[p]>0).mean()*100:.1f}",f"{(mpm.loc[p]>=10).mean()*100:.1f}",f"{mpm.loc[p].median():.2f}",int(rk[p]),int(m.loc[p].max()),int((m.loc[p]>=2).sum()))
        else: w(g,p.replace("hsa-",""),"NA")
S=["hsa-let-7g-5p","hsa-miR-24-3p","hsa-miR-142-3p","hsa-miR-223-3p","hsa-miR-138-5p","hsa-miR-32-3p","hsa-miR-4455","hsa-miR-766-5p","hsa-miR-150-5p","hsa-miR-21-5p","hsa-miR-155-5p"]
for s in m.columns:
    w("SAMPLE", s, int(m[s].sum()), *[f"{np.log2(mpm.loc[p,s]+1):.2f}" if p in m.index else "NA" for p in S])
o.close(); print(open("ta_results.txt").read()[:3000])
