import numpy as np
from scipy.stats import mannwhitneyu
D="/home/claude/ms/data/"
qc={r.split(",")[0]:[float(x) for x in r.split(",")[1:]] for r in open(D+"qc_full.txt").read().strip().split(";")}
comp={r.split(",")[0]:[float(x) for x in r.split(",")[1:]] for r in open(D+"comp.txt").read().strip().split(";")}
S=sorted(qc); y=np.array([s.startswith("BC") for s in S])
g=lambda j:np.array([qc[s][j] for s in S])
reads,kept,noad,badl,ali,umi,det=[g(j) for j in range(7)]
rows=[("Raw reads (millions)",reads/1e6,2),("Reads without adapter (%)",noad/reads*100,1),("Inserts outside 16–30 nt (%)",badl/reads*100,1),
("Retained inserts (millions)",kept/1e6,2),("Retained inserts aligned (%)",ali/kept*100,1),("PCR duplicates (% of aligned)",(1-umi/ali)*100,1),
("UMI molecules (millions)",umi/1e6,2),("miRNAs detected (≥1 molecule)",det,0),
("miR-16-5p (% of molecules)",np.array([comp[s][2] for s in S]),1),("miR-486-5p (% of molecules)",np.array([comp[s][3] for s in S]),1),("miR-122-5p (% of molecules)",np.array([comp[s][4] for s in S]),2)]
out=[]
for lab,v,d in rows:
    f=lambda w:f"{np.median(w):.{d}f} ({np.percentile(w,25):.{d}f}–{np.percentile(w,75):.{d}f})"
    p=mannwhitneyu(v[y],v[~y]).pvalue
    ps="< 0.001" if p<0.001 else f"{p:.3f}" if p<0.1 else f"{p:.2f}"
    out.append([lab,f(v),f(v[y]),f(v[~y]),ps]); print(out[-1])
print("totals: reads",reads.sum()/1e9,"umi",umi.sum()/1e6,"min umi",umi.min(),umi.argmin(),S[umi.argmin()],"max",umi.max(),S[umi.argmax()])
print("libs <500k umi",(umi<5e5).sum(),"p10",np.percentile(umi,10))
import json; json.dump(out,open("qctab.json","w"),ensure_ascii=False)
