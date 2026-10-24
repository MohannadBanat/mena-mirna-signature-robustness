#!/usr/bin/env bash
# Fetch PRJNA1436574 (MENA circulating miRNA, Vishnubalaji et al., J Transl Med 2026)
# macOS (Intel or Apple Silicon). Run in Terminal:  bash fetch_PRJNA1436574.sh [step]
#   step 1 = metadata only (a few KB)      <- run this first, then stop and look
#   step 2 = download all FASTQ (~60-90 GB)
#   step 3 = verify checksums
set -euo pipefail

PROJ="PRJNA1436574"
ROOT="$(cd "$(dirname "$0")" && pwd)/${PROJ}"
META="${ROOT}/metadata"
RAW="${ROOT}/fastq"
STEP="${1:-1}"
mkdir -p "$META" "$RAW"

need() { command -v "$1" >/dev/null 2>&1 || { echo "Installing $2 via Homebrew"; brew install "$2"; }; }

# ---------- STEP 1: metadata ----------
if [[ "$STEP" == "1" ]]; then
  command -v brew >/dev/null || { echo "Homebrew missing: https://brew.sh"; exit 1; }
  need jq jq
  need csvkit csvkit 2>/dev/null || true

  # ENA mirrors SRA and gives direct FASTQ links + MD5s (faster than fasterq-dump)
  FIELDS="run_accession,sample_accession,experiment_accession,sample_alias,sample_title,experiment_title,library_name,instrument_model,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp"
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJ}&result=read_run&fields=${FIELDS}&format=tsv&download=true" \
    -o "${META}/ena_runs.tsv" || echo "ENA not ready yet; falling back to NCBI only"

  # NCBI run table
  curl -fsSL "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=${PROJ}" -o "${META}/ncbi_runinfo.csv" || true

  # BioSample attributes = where clinical labels live, IF the authors submitted them
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJ}&result=read_run&fields=run_accession,sample_accession&format=tsv" \
    | tail -n +2 | cut -f2 | sort -u > "${META}/samples.txt" || true
  : > "${META}/biosample_attributes.jsonl"
  while read -r S; do
    [[ -z "$S" ]] && continue
    curl -fsSL "https://www.ebi.ac.uk/biosamples/samples/${S}.json" \
      | jq -c '{acc: .accession, name: .name, chars: (.characteristics | with_entries(.value = (.value[0].text)))}' \
      >> "${META}/biosample_attributes.jsonl" || echo "fail ${S}"
    sleep 0.2
  done < "${META}/samples.txt"

  echo; echo "=== SUMMARY ==="
  echo "Runs:            $(($(wc -l < "${META}/ena_runs.tsv") - 1))"
  echo "Total FASTQ GB:  $(awk -F'\t' 'NR>1{split($11,a,";"); for(i in a) s+=a[i]} END{printf "%.1f", s/1e9}' "${META}/ena_runs.tsv")"
  echo "Sample titles (first 10 distinct prefixes):"
  awk -F'\t' 'NR>1{print $5}' "${META}/ena_runs.tsv" | sed -E 's/[0-9]+$//' | sort | uniq -c | head
  echo "BioSample attribute keys present:"
  jq -r '.chars | keys[]' "${META}/biosample_attributes.jsonl" | sort | uniq -c
  echo; echo "Free disk: $(df -h "$HOME" | awk 'NR==2{print $4}')"
  echo "Metadata in ${META}. Upload ena_runs.tsv and biosample_attributes.jsonl to Claude before step 2."
  exit 0
fi

# ---------- STEP 2: download FASTQ ----------
if [[ "$STEP" == "2" ]]; then
  need aria2c aria2
  awk -F'\t' 'NR>1{n=split($13,f,";"); for(i=1;i<=n;i++) print "https://"f[i]}' "${META}/ena_runs.tsv" > "${META}/urls.txt"
  echo "Downloading $(wc -l < "${META}/urls.txt") files to ${RAW} (resumable; re-run if interrupted)"
  aria2c -i "${META}/urls.txt" -d "$RAW" -c -x 4 -j 4 --auto-file-renaming=false --summary-interval=60
  exit 0
fi

# ---------- STEP 3: verify ----------
if [[ "$STEP" == "3" ]]; then
  awk -F'\t' 'NR>1{n=split($13,f,";"); split($12,m,";"); for(i=1;i<=n;i++){k=split(f[i],p,"/"); print m[i]"  "p[k]}}' \
    "${META}/ena_runs.tsv" > "${RAW}/md5.txt"
  cd "$RAW"
  bad=0
  while read -r sum file; do
    [[ -f "$file" ]] || { echo "MISSING $file"; bad=1; continue; }
    [[ "$(md5 -q "$file")" == "$sum" ]] || { echo "BAD $file"; bad=1; }
  done < md5.txt
  [[ $bad == 0 ]] && echo "All files verified." || echo "Re-run step 2 for missing/bad files."
fi
