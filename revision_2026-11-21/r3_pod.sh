#!/usr/bin/env bash
# Plan v5, R3 / amendment v5a, A2: alignment sensitivity on a RunPod CPU pod.
# Usage: bash r3_pod.sh setup | test | all | matrix | one RUN SAMPLE
# Trimming, UMI capture and insert filter are copied from pipeline/pod_pipeline.sh.
# --strata is omitted: this bowtie 1.3.1 build exits with an error when it is combined with -k 1,
# where the bowtie manual states it has no effect.
set -euo pipefail
W=/workspace/r3; mkdir -p $W/{fastq,umi,logs,ref,counts/A,counts/B,counts/C}; cd $W
export MAMBA_ROOT_PREFIX=/root/micromamba
export PATH=/root/micromamba/envs/r3/bin:$PATH
JOBS=${JOBS:-14}
GH=https://raw.githubusercontent.com/MohannadBanat/mena-mirna-signature-robustness/revision-2026-11

setup() {
  apt-get update -qq && apt-get install -y -qq curl bzip2 >/dev/null
  [ -x /usr/local/bin/micromamba ] || curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /usr/local bin/micromamba
  [ -d /root/micromamba/envs/r3 ] || micromamba create -y -n r3 -c conda-forge -c bioconda \
    python=3.11 umi_tools=1.1.5 bowtie=1.3.1 bowtie2=2.5.4 samtools=1.20 sra-tools=3.1.1 pigz pandas aria2 >/dev/null
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1436574&result=read_run&fields=run_accession,sample_alias,read_count,fastq_md5,fastq_ftp&format=tsv" > runs.tsv
  echo "runs: $(($(wc -l < runs.tsv)-1))"
  curl -fsSL $GH/rerun_2026-10-24/mkref.py -o mkref.py && python3 mkref.py
  awk '/^>/{keep=($1 ~ /^>hsa-/); if(keep) print $1; next} keep{gsub("U","T"); print}' ref/mature.fa > ref/hsa_mature.fa
  echo "hsa mature: $(grep -c '>' ref/hsa_mature.fa)" | tee ref/n_seq.txt
  sha256sum ref/hsa_mature.fa | tee ref/sha256.txt
  bowtie-build -q ref/hsa_mature.fa ref/bt1 && bowtie2-build -q ref/hsa_mature.fa ref/bt2
  # versions
  { bowtie --version | head -1; bowtie2 --version | head -1; samtools --version | head -1; umi_tools --version; } > logs/versions.txt 2>&1
  cat logs/versions.txt; df -h /workspace | tail -1
}

one() {  # $1 = run accession, $2 = sample alias
  local R=$1 S=$2
  [ -s counts/C/$S.umi.tsv ] && return 0
  local URL MD5
  URL=$(awk -F'\t' -v r=$R '$1==r{print "https://"$5}' runs.tsv)
  MD5=$(awk -F'\t' -v r=$R '$1==r{print $4}' runs.tsv)
  # ENA served this pod at ~0.15 MB/s on 22 Nov 2026, so reads come from NCBI's SRA mirror on AWS instead;
  # completeness is checked afterwards against ENA's read_count (runs.tsv), since the ENA MD5 does not apply.
  if [ ! -s fastq/$R.fastq.gz ]; then
    curl -fsSL --retry 5 -o fastq/$R.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/$R/$R
    fasterq-dump -e 2 --split-spot --stdout fastq/$R.sra 2> logs/$S.fqd.txt | pigz -1 -p 2 > fastq/$R.fastq.gz
    rm -f fastq/$R.sra
  fi
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
            m=re.search("AACTGTAGGCAC|ACTGTAGGCACC",s)
            if not m: noad+=1; continue
            i=m.start()-(1 if m.group(0).startswith("ACTG") else 0)
        ins=s[:i]; umi=s[i+len(ad):i+len(ad)+12]
        if len(umi)<12 or not (16<=len(ins)<=30): short+=1; continue
        o.write(f"@{h[1:].split()[0]}_{umi}\n{ins}\n+\n{q[:i]}\n"); kept+=1
with open(f"logs/{S}.trim.txt","w") as L:
    L.write(f"{S}\t{R}\t{n}\t{kept}\t{noad}\t{short}\n")
PY
  rm -f fastq/$R.fastq.gz
  for X in A B C; do
    case $X in
      A) bowtie -p 1 -n 1 -l 16 --best -k 1 --norc -S ref/bt1 umi/$S.fq 2> logs/$S.A.txt ;;
      B) bowtie -p 1 -v 0 --best -k 1 --norc -S ref/bt1 umi/$S.fq 2> logs/$S.B.txt ;;
      C) bowtie2 -p 1 --end-to-end --very-sensitive -L 10 -N 0 --norc --seed 0 -x ref/bt2 -U umi/$S.fq 2> logs/$S.C.txt ;;
    esac | samtools view -b -F 4 - | samtools sort -o umi/$S.$X.bam - && samtools index umi/$S.$X.bam
    samtools idxstats umi/$S.$X.bam | awk '$3>0{print $1"\t"$3}' > counts/$X/$S.reads.tsv
    umi_tools dedup -I umi/$S.$X.bam --per-contig --method=unique --umi-separator=_ -S umi/$S.$X.dedup.bam -L logs/$S.$X.dedup.txt
    samtools index umi/$S.$X.dedup.bam
    samtools idxstats umi/$S.$X.dedup.bam | awk '$3>0{print $1"\t"$3}' > counts/$X/$S.umi.tsv.tmp
    mv counts/$X/$S.umi.tsv.tmp counts/$X/$S.umi.tsv
    rm -f umi/$S.$X.bam* umi/$S.$X.dedup.bam*
  done
  rm -f umi/$S.fq
}

case "${1:-}" in
  setup) setup ;;
  one) one "$2" "$3" ;;
  test)  R=$(awk -F'\t' 'NR==2{print $1}' runs.tsv); S=$(awk -F'\t' 'NR==2{print $2}' runs.tsv)
         time one $R $S; cat logs/$S.trim.txt
         for X in A B C; do echo "$X: $(wc -l < counts/$X/$S.umi.tsv) miRNAs, $(awk '{s+=$2}END{print s}' counts/$X/$S.umi.tsv) molecules"; done ;;
  all)   awk -F'\t' 'NR>1{print $1" "$2}' runs.tsv | xargs -P $JOBS -n 2 bash $0 one || true
         echo "done: $(ls counts/C/*.umi.tsv | wc -l) of 200"; echo ALLDONE ;;
  matrix) python3 - <<'PY'
import pandas as pd,glob,os
for X in "ABC":
    for u in ["umi","reads"]:
        d={os.path.basename(f).split(".")[0]:pd.read_csv(f,sep="\t",header=None,index_col=0)[1] for f in glob.glob(f"counts/{X}/*.{u}.tsv")}
        m=pd.DataFrame(d).fillna(0).astype(int).sort_index(axis=1); m.to_csv(f"r3_{X}_{u}.tsv",sep="\t"); print(X,u,m.shape,int(m.values.sum()))
q=pd.concat([pd.read_csv(f,sep="\t",header=None) for f in glob.glob("logs/*.trim.txt")]); q.columns=["sample","run","reads","kept","no_adapter","bad_len"]
q.sort_values("sample").to_csv("r3_qc.tsv",sep="\t",index=False)
PY
         tar czf r3_results.tar.gz r3_*.tsv ref/n_seq.txt ref/sha256.txt logs/versions.txt
         sha256sum r3_results.tar.gz; ls -la r3_results.tar.gz ;;
  *) echo "usage: $0 setup|test|all|matrix|one R S" ;;
esac
