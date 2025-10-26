#!/usr/bin/env bash
# ===============================================================
# Script Name: 02_ldsc_heritability_qc.sh
# Description: Perform heritability QC (h²ₛₙₚ > 0.05) for base data using LDSC.
# Usage: bash 02_ldsc_heritability_qc.sh
# ===============================================================

# QC control for base data
# 1. Heritability check: h²ₛₙₚ > 0.05

# Executable in any directory
module load python/3.9.7 2>/dev/null || module load python 2>/dev/null || true
module load git 2>/dev/null || true

# Create and activate virtual environment (under $WORK)
export LDSC_ENV="$WORK/ldsc_env"
python3 -m venv "$LDSC_ENV"
source "$LDSC_ENV/bin/activate"

# Install dependencies
pip install -U pip
pip install "numpy==1.26.4" "pandas==2.1.4" "scipy==1.11.4" "numexpr==2.8.7" ldsc

# Self-check (help message indicates success)
which munge_sumstats.py && munge_sumstats.py --help | head -n 3
which ldsc.py           && ldsc.py           --help | head -n 3

# Reference directory
export REF_DIR="${WORK:-$HOME}/ldsc_ref"
cd "$REF_DIR"

# Upload local reference files to TACC
scp -r /Users/li/Desktop/UTAustin/Lab_Rotation/Arbel_Lab/Project/reference/eur_w_ld_chr.tar.gz liliu@ls6.tacc.utexas.edu:/work/11063/liliu/ls6/ldsc_ref
scp -r /Users/li/Desktop/UTAustin/Lab_Rotation/Arbel_Lab/Project/reference/w_hm3.snplist.gz liliu@ls6.tacc.utexas.edu:/work/11063/liliu/ls6/ldsc_ref

# Decompress reference data
gunzip -f w_hm3.snplist.gz 
tar -xzf eur_w_ld_chr.tar.gz

# Return to base data directory

# Generate .sumstats.gz
export REF_DIR="${WORK:-$HOME}/ldsc_ref"

for f in *_base_prsice.txt; do
  p="${f%_base_prsice.txt}"
  echo "==> $p"
  munge_sumstats.py \
    --sumstats "$f" \
    --snp SNP --a1 A1 --a2 A2 \
    --p P --signed-sumstats BETA,0 \
    --N 22497 \
    --no-alleles \
    --ignore INFO,MAF \
    --info-min 0 --maf-min 0 \
    --chunksize 200000 \
    --out "$p" || { echo "[munge failed] $p"; continue; }

  ldsc.py \
    --h2 "${p}.sumstats.gz" \
    --ref-ld-chr "$REF_DIR/eur_w_ld_chr/" \
    --w-ld-chr  "$REF_DIR/eur_w_ld_chr/" \
    --out "${p}.ldsc" || echo "[ldsc failed] $p"
done

# Display diagnostics
for log in *.ldsc.log; do
  [ -e "$log" ] || { echo "No *.ldsc.log available"; break; }
  echo "=== ${log%.ldsc.log} ==="
  grep -E "Total Observed scale h2:|Mean chi|Intercept|Ratio" "$log" | sed 's/^/  /' || echo "  [no diagnostics]"
  echo
done

# Summarize results (h2, SE, Intercept, MeanChi2, Ratio)
out=ldsc_h2_summary.tsv
echo -e "GWAS\tObserved_h2\tSE\tIntercept\tMeanChi2\tRatio" > "$out"

awk '
  BEGIN{ OFS="\t" }
  FILENAME ~ /\.ldsc\.log$/ {
    if ($0 ~ /Total Observed scale h2:/ && match($0,/h2:\s*([-+]?[0-9.]+)\s*\(([-+]?[0-9.]+)\)/,a)) {h2=a[1]; se=a[2]}
    if ($0 ~ /^Intercept:/                   && match($0,/Intercept:\s*([-+]?[0-9.]+)/,b)) {ic=b[1]}
    if ($0 ~ /^Mean chi\^?2:/                && match($0,/: *([-+]?[0-9.]+)/,c))           {mc=c[1]}
    if ($0 ~ /^Ratio/                        && match($0,/Ratio:\s*([-+]?[0-9.]+)/,d))     {ra=d[1]}
  }
  ENDFILE {
    gw=FILENAME; sub(/\.ldsc\.log$/, "", gw)
    if (h2=="") h2="NA"; if (se=="") se="NA"; if (ic=="") ic="NA"; if (mc=="") mc="NA"; if (ra=="") ra="NA";
    print gw, h2, se, ic, mc, ra;
    h2=se=ic=mc=ra="";
  }
' ieu-b-48*.ldsc.log >> "$out"

column -t ldsc_h2_summary.tsv | less -S

# Filter h² > 0.05 and sort descending
awk -F'\t' 'NR==1 || ($2+0)>0.05' ldsc_h2_summary.tsv > ldsc_h2_gt0.05.tsv
column -t ldsc_h2_gt0.05.tsv | less -S
{ head -n1 ldsc_h2_gt0.05.tsv; tail -n +2 ldsc_h2_gt0.05.tsv | sort -k2,2nr; } | column -t

# Organize output files
outdir="base_gwas_h2_0.05"
mkdir -p "$outdir"
awk -F'\t' 'NR>1 {print $1}' ldsc_h2_gt0.05.tsv > keep_ids.txt

while read -r id; do
  for f in \
    "${id}_base_prsice.txt" \
    "${id}.sumstats.gz" \
    "${id}.sumstats.gz.tbi" \
    "${id}.ldsc.log" \
    "${id}.vcf.gz"
  do
    [ -f "$f" ] && cp -p "$f" "$outdir"/
  done
done < keep_ids.txt

ls -lh "$outdir"
