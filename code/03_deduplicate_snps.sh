#!/usr/bin/env bash
# ===============================================================
# Script Name: 03_deduplicate_snps.sh
# Description: Remove duplicate SNPs and keep entries with the best (lowest) P-value.
# Usage: bash 03_deduplicate_snps.sh
# ===============================================================

set -euo pipefail

outdir="base_dedup_bestP"
mkdir -p "$outdir"

summary="dedup_summary.tsv"
echo -e "GWAS\tInput_rows\tRemoved_SNPdot\tRemoved_duplicates\tKept_rows\tTotal_removed" > "$summary"

for in in *_base_prsice.txt; do
  [ -f "$in" ] || { echo "No *_base_prsice.txt files found"; break; }

  id="${in%_base_prsice.txt}"
  clean="${outdir}/${id}_base_prsice.clean.txt"
  out="${outdir}/${id}_base_prsice.bestP.nodup.txt"

  echo "==> $id  (drop SNP='.' → dedup best P)"

  in_rows=$(( $(wc -l < "$in") - 1 ))

  awk -F'\t' 'NR==1{print; next} $3!="."' "$in" > "$clean"
  clean_rows=$(( $(wc -l < "$clean") - 1 ))
  removed_dot=$(( in_rows - clean_rows ))

  awk -F'\t' '
    NR==1{hdr=$0; next}
    {
      key=$3
      p  = ($8==""||$8=="."||$8=="NA") ? 1e300 : $8+0
      se = ($7==""||$7=="."||$7=="NA") ? 1e300 : $7+0
      if(!(key in best) || p < bestp[key] || (p==bestp[key] && se < bestse[key])) {
        best[key]=$0; bestp[key]=p; bestse[key]=se
      }
    }
    END{
      print hdr
      for(k in best) print best[k]
    }
  ' "$clean" | ( read -r h; echo "$h"; sort -k1,1 -k2,2n ) > "$out"

  kept_rows=$(( $(wc -l < "$out") - 1 ))
  removed_dup=$(( clean_rows - kept_rows ))
  total_removed=$(( in_rows - kept_rows ))

  printf "%s\t%d\t%d\t%d\t%d\t%d\n" \
    "$id" "$in_rows" "$removed_dot" "$removed_dup" "$kept_rows" "$total_removed" >> "$summary"
done
