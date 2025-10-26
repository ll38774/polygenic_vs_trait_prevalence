#!/usr/bin/env bash
# ===============================================================
# Script Name: 04_remove_ambiguous_snps.sh
# Description: Remove non-ACGT and palindromic SNPs (A/T, T/A, C/G, G/C).
# Usage: bash 04_remove_ambiguous_snps.sh
# ===============================================================

# 1) Enter current directory containing *_base_prsice.bestP.nodup.txt
# 2) Create output directory and summary table
mkdir -p base_qc_nopal
echo -e "GWAS\tInput_rows\tRemoved_nonACGT\tRemoved_palindromic\tKept_rows" > base_qc_nopal/nopal_summary.tsv

# 3) Batch process all *_base_prsice.bestP.nodup.txt files
for in in *_base_prsice.bestP.nodup.txt; do
  [ -f "$in" ] || { echo "No *_base_prsice.bestP.nodup.txt found"; break; }
  id="${in%_base_prsice.bestP.nodup.txt}"
  out="base_qc_nopal/${id}_base_prsice.bestP.nodup.nopal.txt.gz"

  inp=$(( $(wc -l < "$in") - 1 ))

  { 
    head -n1 "$in"
    awk -F'\t' '
      $4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ && !($4$5~/^(AT|TA|CG|GC)$/)
    ' "$in" | sed -n '2,$p'
  } | gzip > "$out"

  nonacgt=$(awk -F'\t' 'NR>1 && !($4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/){c++} END{print c+0}' "$in")
  pal=$(awk -F'\t' 'NR>1 && ($4$5=="AT"||$4$5=="TA"||$4$5=="CG"||$4$5=="GC"){c++} END{print c+0}' "$in")
  kept=$(( inp - nonacgt - pal ))

  printf "%s\t%d\t%d\t%d\t%d\n" "$id" "$inp" "$nonacgt" "$pal" "$kept" >> base_qc_nopal/nopal_summary.tsv
  echo "nopal -> $out   (kept ${kept}/${inp}, removed $((nonacgt+pal)))"
done

# 4) Quick preview
column -t base_qc_nopal/nopal_summary.tsv | head
ls -lh base_qc_nopal/*nopal.txt.gz | head
