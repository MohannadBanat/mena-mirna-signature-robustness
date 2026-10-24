import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
S=[x.split(',') for x in open('data/tissue_samples.txt').read().strip().split(';')]
cols=["let-7g-5p","miR-24-3p","miR-142-3p","miR-223-3p","miR-138-5p","miR-32-3p","miR-4455","miR-766-5p","miR-150-5p","miR-21-5p","miR-155-5p"]
X={c:np.array([float(r[i+2]) for r in S]) for i,c in enumerate(cols)}
for a in ["miR-142-3p","miR-223-3p"]:
    for b in ["miR-150-5p","miR-155-5p","miR-21-5p"]:
        r,p=spearmanr(X[a],X[b]); print(f"Q3 {a} vs {b}: rho={r:.3f} p={p:.2g}")
r,p=spearmanr(X["miR-142-3p"],X["miR-223-3p"]); print(f"142 vs 223 rho={r:.3f} p={p:.2g}")
T={}
for l in open('data/tissue_summary.tsv'):
    p=l.rstrip('\n').split('\t')
    if p[0] in ('UP','DN','PANEL','REF') and len(p)>3: T[(p[0],p[1])]=[float(x) for x in p[2:]]
up=[v[2] for (g,m),v in T.items() if g=='UP']; dn=[v[2] for (g,m),v in T.items() if g=='DN']
print("Q4 up n",len(up),"median MPM",np.median(up),"IQR",np.percentile(up,[25,75]),"| down n",len(dn),"median",np.median(dn),np.percentile(dn,[25,75]),"MWU p",mannwhitneyu(up,dn).pvalue)
upd=[v[0] for (g,m),v in T.items() if g=='UP']; dnd=[v[0] for (g,m),v in T.items() if g=='DN']
print("detect>=1 up median",np.median(upd),"rare(<50%) up:",sum(d<50 for d in upd),"of",len(upd),"; down rare:",sum(d<50 for d in dnd),"of",len(dnd))
print("rare up list",[m for (g,m),v in T.items() if g=='UP' and v[0]<50])
# Q1 classification
for (g,m),v in T.items():
    if g=='PANEL':
        cls="supported" if v[1]>=50 else ("rare" if v[0]<50 else "low")
        print("Q1",m,v[0],v[1],v[2],int(v[3]),cls)
# plasma DE miRNAs: tissue abundance vs plasma log2FC (exploratory)
d={r.split(',')[0].replace('hsa-',''):r.split(',')[1:] for r in open('data/de.txt').read().strip().split(';')}
pairs=[(float(d[m][1]),v[2]) for (g,m),v in T.items() if g in('UP','DN') and m in d]
r,p=spearmanr([a for a,b in pairs],[np.log2(b+1) for a,b in pairs]); print("expl: plasma lfc vs tissue log MPM among 37 DE: rho",round(r,3),"p",p)
print("4455 samples detected:",sum(X['miR-4455']>0), "766:",sum(X['miR-766-5p']>0))
