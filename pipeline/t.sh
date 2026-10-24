#!/usr/bin/env bash
# PRJNA953015 (Elango 2023, 96 FFPE tumours, QIAseq miRNA) - same processing as the MENA plasma pipeline
set -euo pipefail
W=/root/tissue; mkdir -p $W/{fastq,umi,counts,logs,ref}; cd $W
PROJ=PRJNA953015; THREADS=${THREADS:-2}
export MAMBA_ROOT_PREFIX=/root/micromamba; export PATH=/root/micromamba/envs/mirna/bin:$PATH
setup() {
  apt-get update -qq && apt-get install -y -qq curl bzip2 >/dev/null
  [ -x /usr/local/bin/micromamba ] || curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /usr/local bin/micromamba
  [ -d /root/micromamba/envs/mirna ] || micromamba create -y -n mirna -c conda-forge -c bioconda python=3.11 umi_tools=1.1.5 bowtie=1.3.1 samtools=1.20 pandas aria2 >/dev/null
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$PROJ&result=read_run&fields=run_accession,sample_title,read_count,fastq_md5,fastq_ftp&format=tsv" > runs.tsv
  curl -fsSL https://www.mirbase.org/download/mature.fa -o ref/mature.fa || curl -fsSL https://mirbase.org/ftp/22.1/mature.fa.gz | gunzip > ref/mature.fa
  curl -fsSL https://mirbase.org/ftp/22.1/mature.fa.gz | gunzip > ref/mature22.fa || true
  awk '/^>/{keep=($1 ~ /^>hsa-/); if(keep) print $1; next} keep{gsub("U","T"); print}' ref/mature.fa > ref/hsa_mature.fa
  bowtie-build -q ref/hsa_mature.fa ref/hsa_mature
  echo "runs: $(($(wc -l < runs.tsv)-1))  hsa: $(grep -c '>' ref/hsa_mature.fa)  v22: $(grep -c '>hsa' ref/mature22.fa)"; df -h /root | tail -1
}
one() {
  local R=$1 S=$2; [ -s counts/$S.tsv ] && return 0
  local URL MD5; URL=$(awk -F'\t' -v r=$R '$1==r{print "https://"$5}' runs.tsv); MD5=$(awk -F'\t' -v r=$R '$1==r{print $4}' runs.tsv)
  [ -s fastq/$R.fastq.gz ] || aria2c -q -x 4 -d fastq -o $R.fastq.gz "$URL"
  echo "$MD5  fastq/$R.fastq.gz" | md5sum -c --quiet - || { echo "MD5 FAIL $R"; rm -f fastq/$R.fastq.gz; return 1; }
  python3 - "$R" "$S" <<'PY'
import gzip,sys,re
R,S=sys.argv[1],sys.argv[2]; ad="AACTGTAGGCACCATCAAT"; n=kept=noad=short=0
with gzip.open(f"fastq/{R}.fastq.gz","rt") as f, open(f"umi/{S}.fq","w") as o:
    while True:
        h=f.readline()
        if not h: break
        s=f.readline().strip(); f.readline(); q=f.readline().strip(); n+=1
        i=s.find(ad)
        if i<0:
            m=re.search("AACTGTAGGCAC|ACTGTAGGCACC",s)
            if not m: noad+=1; continue
            i=m.start()-(1 if m.group(0).startswith("ACTG") else 0)
        ins=s[:i]; umi=s[i+len(ad):i+len(ad)+12]
        if len(umi)<12 or not (16<=len(ins)<=30): short+=1; continue
        o.write(f"@{h[1:].split()[0]}_{umi}\n{ins}\n+\n{q[:i]}\n"); kept+=1
open(f"logs/{S}.trim.txt","w").write(f"{S}\t{R}\t{n}\t{kept}\t{noad}\t{short}\n")
PY
  bowtie -p $THREADS -v 1 --best -k 1 --norc -S ref/hsa_mature umi/$S.fq 2> logs/$S.bowtie.txt | samtools view -b -F 4 - | samtools sort -o umi/$S.bam - && samtools index umi/$S.bam
  umi_tools dedup -I umi/$S.bam --per-gene --per-contig --method=unique --umi-separator=_ -S umi/$S.dedup.bam -L logs/$S.dedup.txt
  samtools idxstats umi/$S.bam | awk '$3>0{print $1"\t"$3}' > counts/$S.raw.tsv
  samtools index umi/$S.dedup.bam; samtools idxstats umi/$S.dedup.bam | awk '$3>0{print $1"\t"$3}' > counts/$S.tsv
  rm -f umi/$S.fq umi/$S.bam* umi/$S.dedup.bam* fastq/$R.fastq.gz
}
export -f one; export W THREADS
case "${1:-}" in
  setup) setup ;;
  test) R=$(awk -F'\t' 'NR==2{print $1}' runs.tsv); S=$(awk -F'\t' 'NR==2{print $2}' runs.tsv); time one $R $S; cat logs/$S.trim.txt; tail -2 logs/$S.bowtie.txt; echo "molecules: $(awk '{s+=$2}END{print s}' counts/$S.tsv) miRNAs: $(wc -l < counts/$S.tsv)"; sort -k2,2nr counts/$S.tsv | head -8 ;;
  all) awk -F'\t' 'NR>1{print $1" "$2}' runs.tsv | xargs -P ${JOBS:-4} -n 2 bash -c 'cd $W; one "$0" "$1" || echo "FAILED $1" >> logs/failed.txt'; echo "DONE $(ls counts/*.tsv | grep -vc raw)" ;;
  matrix) python3 - <<'PY'
import pandas as pd,glob,os
mat={}
for l in open("ref/mature.fa"):
    if l.startswith(">hsa-"): p=l[1:].split(); mat[p[0]]=p[1]
v22={}
for l in open("ref/mature22.fa"):
    if l.startswith(">hsa-"): p=l[1:].split(); v22[p[1]]=p[0]
d={os.path.basename(f)[:-4]:pd.read_csv(f,sep="\t",header=None,index_col=0)[1] for f in glob.glob("counts/*.tsv") if not f.endswith(".raw.tsv")}
m=pd.DataFrame(d).fillna(0).astype(int).sort_index(axis=1)
m.index=[v22.get(mat.get(i,""),i) for i in m.index]; m=m.groupby(level=0).sum()
m.to_csv("tissue_counts_umi_v22.tsv",sep="\t")
q=pd.concat([pd.read_csv(f,sep="\t",header=None) for f in glob.glob("logs/*.trim.txt")]); q.columns=["sample","run","reads","kept","no_adapter","bad_len"]
q["umi_total"]=q["sample"].map(m.sum()); q["mirnas_detected"]=q["sample"].map((m>0).sum()); q.sort_values("sample").to_csv("tissue_qc.tsv",sep="\t",index=False)
mpm=m.div(m.sum(),axis=1)*1e6
P=["hsa-let-7g-5p","hsa-miR-138-5p","hsa-miR-142-3p","hsa-miR-24-3p","hsa-miR-32-3p","hsa-miR-223-3p","hsa-miR-4455","hsa-miR-766-5p","hsa-miR-106b-3p","hsa-miR-4497","hsa-miR-660-3p","hsa-miR-30a-5p","hsa-miR-3195","hsa-miR-500b-5p","hsa-miR-1976","hsa-miR-106a-3p","hsa-miR-382-3p","hsa-miR-4639-5p","hsa-miR-6741-3p","hsa-miR-150-5p","hsa-miR-23a-3p","hsa-miR-451a","hsa-miR-16-5p","hsa-miR-122-5p","hsa-miR-486-5p"]
rk=mpm.median(1).rank(ascending=False)
with open("tissue_panel_summary.txt","w") as o:
    o.write(f"matrix {m.shape}; molecules median {int(m.sum().median())}; miRNAs/library median {int((m>0).sum().median())}\n")
    o.write("mirna\tpct_ge1mol\tpct_ge10mpm\tmedian_mpm\trank\n")
    for p in P:
        if p in m.index: o.write(f"{p}\t{(m.loc[p]>0).mean()*100:.1f}\t{(mpm.loc[p]>=10).mean()*100:.1f}\t{mpm.loc[p].median():.2f}\t{int(rk[p])}\n")
        else: o.write(f"{p}\tNA\tNA\tNA\tNA\n")
print(open("tissue_panel_summary.txt").read()); print(q.describe().round(0))
PY
  tar czf tissue_results.tar.gz tissue_counts_umi_v22.tsv tissue_qc.tsv tissue_panel_summary.txt runs.tsv; ls -la tissue_results.tar.gz ;;
esac
