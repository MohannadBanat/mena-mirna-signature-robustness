#!/usr/bin/env bash
# PRJNA1436574 QIAseq miRNA pipeline for a RunPod Ubuntu CPU pod.
# Usage:  bash pod_pipeline.sh setup | test | all | matrix
#   setup  : install tools (micromamba env), fetch run list + miRBase mature human
#   test   : run ONE sample end-to-end and print QC  <- check before 'all'
#   all    : download + process all 200 samples (resumable; skips finished ones)
#   matrix : merge per-sample counts into counts_matrix.tsv + qc_table.tsv
set -euo pipefail
W=/root/mena; mkdir -p $W/{fastq,trim,umi,counts,logs,ref}
cd $W
PROJ=PRJNA1436574
THREADS=${THREADS:-8}
ENVBIN=/root/micromamba/envs/mirna/bin
export PATH=$ENVBIN:$PATH

# QIAseq miRNA: read = miRNA + 3' adapter AACTGTAGGCACCATCAAT + 12-nt UMI (+ more adapter)
ADAPTER=AACTGTAGGCACCATCAAT

setup() {
  apt-get update -qq && apt-get install -y -qq curl bzip2 pigz >/dev/null
  if [ ! -x /usr/local/bin/micromamba ]; then
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /usr/local bin/micromamba
  fi
  export MAMBA_ROOT_PREFIX=/root/micromamba
  micromamba create -y -n mirna -c conda-forge -c bioconda \
    python=3.11 cutadapt=4.9 umi_tools=1.1.5 bowtie=1.3.1 samtools=1.20 pandas aria2 >/dev/null
  # run list with ftp paths + md5
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$PROJ&result=read_run&fields=run_accession,sample_alias,read_count,fastq_md5,fastq_ftp&format=tsv" > runs.tsv
  echo "runs: $(($(wc -l < runs.tsv)-1))"
  # miRBase 22.1 mature, human only, U->T
  curl -fsSL https://www.mirbase.org/download/mature.fa -o ref/mature.fa || \
    curl -fsSL https://mirbase.org/ftp/22.1/mature.fa.gz | gunzip > ref/mature.fa
  awk '/^>/{keep=($1 ~ /^>hsa-/); if(keep) print $1; next} keep{gsub("U","T"); print}' ref/mature.fa > ref/hsa_mature.fa
  echo "hsa mature: $(grep -c '>' ref/hsa_mature.fa)"
  bowtie-build -q ref/hsa_mature.fa ref/hsa_mature
  df -h /root | tail -1
}

one() {  # $1 = run accession, $2 = sample alias
  local R=$1 S=$2
  [ -s counts/$S.tsv ] && return 0
  local URL MD5
  URL=$(awk -F'\t' -v r=$R '$1==r{print "https://"$5}' runs.tsv)
  MD5=$(awk -F'\t' -v r=$R '$1==r{print $4}' runs.tsv)
  [ -s fastq/$R.fastq.gz ] || aria2c -q -x 4 -d fastq -o $R.fastq.gz "$URL"
  echo "$MD5  fastq/$R.fastq.gz" | md5sum -c --quiet - || { echo "MD5 FAIL $R"; rm -f fastq/$R.fastq.gz; return 1; }
  # 1) cut adapter, keep the 12 nt that follow it as UMI -> put UMI into read name
  python3 - "$R" "$S" <<'PY'
import gzip,sys,re
R,S=sys.argv[1],sys.argv[2]
ad="AACTGTAGGCACCATCAAT"
n=kept=noad=short=0
with gzip.open(f"fastq/{R}.fastq.gz","rt") as f, open(f"umi/{S}.fq","w") as o:
    while True:
        h=f.readline()
        if not h: break
        s=f.readline().strip(); f.readline(); q=f.readline().strip(); n+=1
        i=s.find(ad)
        if i<0:
            # allow 1 mismatch in first 12 nt of adapter
            m=re.search("AACTGTAGGCAC|ACTGTAGGCACC",s)
            if not m: noad+=1; continue
            i=m.start()-(1 if m.group(0).startswith("ACTG") else 0)
        ins=s[:i]; umi=s[i+len(ad):i+len(ad)+12]
        if len(umi)<12 or not (16<=len(ins)<=30): short+=1; continue
        o.write(f"@{h[1:].split()[0]}_{umi}\n{ins}\n+\n{q[:i]}\n"); kept+=1
with open(f"logs/{S}.trim.txt","w") as L:
    L.write(f"{S}\t{R}\t{n}\t{kept}\t{noad}\t{short}\n")
PY
  # 2) align to hsa mature (end-to-end, <=1 mismatch, best hit)
  bowtie -p $THREADS -v 1 --best -k 1 --norc -S ref/hsa_mature umi/$S.fq 2> logs/$S.bowtie.txt \
    | samtools view -b -F 4 - | samtools sort -o umi/$S.bam - && samtools index umi/$S.bam
  # 3) UMI dedup per miRNA (reference name = feature)
  umi_tools dedup -I umi/$S.bam --per-gene --per-contig --method=unique --umi-separator=_ \
    -S umi/$S.dedup.bam -L logs/$S.dedup.txt
  samtools idxstats umi/$S.bam       | awk '$3>0{print $1"\t"$3}' > counts/$S.raw.tsv
  samtools index umi/$S.dedup.bam
  samtools idxstats umi/$S.dedup.bam | awk '$3>0{print $1"\t"$3}' > counts/$S.tsv
  # 4) free disk
  rm -f umi/$S.fq umi/$S.bam* umi/$S.dedup.bam* fastq/$R.fastq.gz
}

case "${1:-}" in
  setup) setup ;;
  test)  R=$(awk -F'\t' 'NR==2{print $1}' runs.tsv); S=$(awk -F'\t' 'NR==2{print $2}' runs.tsv)
         time one $R $S; cat logs/$S.trim.txt; tail -3 logs/$S.bowtie.txt
         echo "miRNAs detected: $(wc -l < counts/$S.tsv)  UMI-collapsed total: $(awk '{s+=$2}END{print s}' counts/$S.tsv)"
         sort -k2,2nr counts/$S.tsv | head -15 ;;
  all)   awk -F'\t' 'NR>1{print $1"\t"$2}' runs.tsv | while IFS=$'\t' read R S; do
           echo "[$(date +%T)] $S"; one $R $S || echo "FAILED $S" >> logs/failed.txt; done
         echo DONE ;;
  matrix) python3 - <<'PY'
import pandas as pd,glob,os
d={os.path.basename(f)[:-4]:pd.read_csv(f,sep="\t",header=None,index_col=0)[1] for f in glob.glob("counts/*.tsv") if not f.endswith(".raw.tsv")}
r={os.path.basename(f)[:-8]:pd.read_csv(f,sep="\t",header=None,index_col=0)[1] for f in glob.glob("counts/*.raw.tsv")}
m=pd.DataFrame(d).fillna(0).astype(int).sort_index(axis=1); m.to_csv("counts_matrix_umi.tsv",sep="\t")
pd.DataFrame(r).fillna(0).astype(int).sort_index(axis=1).to_csv("counts_matrix_reads.tsv",sep="\t")
q=pd.concat([pd.read_csv(f,sep="\t",header=None) for f in glob.glob("logs/*.trim.txt")])
q.columns=["sample","run","reads","kept","no_adapter","bad_len"]
q["umi_total"]=q["sample"].map(m.sum()); q["mirnas_detected"]=q["sample"].map((m>0).sum())
q.sort_values("sample").to_csv("qc_table.tsv",sep="\t",index=False)
print(m.shape); print(q.describe().round(1))
PY
         tar czf mena_counts.tar.gz counts_matrix_umi.tsv counts_matrix_reads.tsv qc_table.tsv logs/*.trim.txt
         ls -la mena_counts.tar.gz ;;
  *) echo "usage: $0 setup|test|all|matrix" ;;
esac
