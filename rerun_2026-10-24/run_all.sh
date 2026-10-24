#!/usr/bin/env bash
# Full rerun: pipeline (parallel) -> matrix -> analysis.py -> gse.py -> extra.py. Progress in /root/mena/out (served on :8888).
set -uo pipefail
W=/root/mena; mkdir -p $W/out $W/results; cd $W
exec > >(tee -a $W/out/run_all.log) 2>&1
echo "START $(date -u +%FT%TZ)"
( sleep 23400; echo "AUTO-STOP $(date -u +%FT%TZ)" >> $W/out/run_all.log; runpodctl stop pod "$RUNPOD_POD_ID" ) >/dev/null 2>&1 &
( cd $W/out && nohup python3 -m http.server 8888 >/dev/null 2>&1 & )
bash pod_pipeline.sh setup || { echo "SETUP FAILED"; exit 1; }
export PATH=/root/micromamba/envs/mirna/bin:$PATH
pip install -q pydeseq2==0.5.4 scikit-learn scipy 2>&1 | tail -2
echo "SETUP DONE $(date -u +%FT%TZ)"
export THREADS=1
( while sleep 120; do echo "$(date -u +%T) done $(ls counts/*.raw.tsv 2>/dev/null | wc -l)/200 failed $(cat logs/failed.txt 2>/dev/null | wc -l) disk $(df -h /root | awk 'NR==2{print $4}')" > out/progress.txt; done ) &
PROG=$!
awk -F'\t' 'NR>1{print $1" "$2}' runs.tsv | xargs -P 14 -n 2 bash -c 'bash /root/mena/pod_pipeline.sh one "$0" "$1" >/dev/null 2>>/root/mena/logs/one_err.txt || echo "FAILED $1" >> /root/mena/logs/failed.txt'
# one retry pass for failures
if [ -s logs/failed.txt ]; then mv logs/failed.txt logs/failed_pass1.txt
  awk -F'\t' 'NR>1{print $1" "$2}' runs.tsv | xargs -P 6 -n 2 bash -c 'bash /root/mena/pod_pipeline.sh one "$0" "$1" >/dev/null 2>>/root/mena/logs/one_err.txt || echo "FAILED $1" >> /root/mena/logs/failed.txt'; fi
kill $PROG 2>/dev/null
echo "PIPELINE DONE $(date -u +%FT%TZ) samples=$(ls counts/*.raw.tsv | wc -l)"; cat logs/failed.txt 2>/dev/null
bash pod_pipeline.sh matrix | tail -15
N_ITER=2 python analysis.py $W > /dev/null 2>&1; echo "analysis.py exit $?"; cp results/summary.txt out/analysis_summary.txt
python gse.py $W > /dev/null 2>&1; echo "gse.py exit $?"; cp results/gse73002_summary.txt out/ 2>/dev/null
python extra.py $W 2>&1 | tail -3; echo "extra.py exit"; cp results/extra_summary.txt out/ 2>/dev/null
cp results/*.tsv out/ 2>/dev/null; cp counts_matrix_umi.tsv counts_matrix_reads.tsv qc_table.tsv out/ 2>/dev/null
tar czf out/mena_rerun_results.tar.gz results counts_matrix_umi.tsv counts_matrix_reads.tsv qc_table.tsv logs/*.trim.txt 2>/dev/null
echo "ALLDONE $(date -u +%FT%TZ)"
