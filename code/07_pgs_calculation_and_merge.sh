#!/usr/bin/env bash
# ===============================================================
# Script Name: 07_pgs_calculation_and_merge.sh
# Description: Calculate polygenic scores (PGS) for multiple traits across
#              1000 Genomes populations using PLINK 1.9, then merge results
#              with population and trait annotations.
# ===============================================================


# 4) PGS calculation (PLINK)
#!/usr/bin/env bash
set -euo pipefail

PLINK=/work/11063/liliu/ls6/tools/plink1.9/plink
BASE_DIR=/work/11063/liliu/ls6/clean_data/opengwas_clean_data
TGT_DIR=/work/11063/liliu/ls6/clean_data/1000_genome_clean_data_final
OUTDIR=/work/11063/liliu/ls6/pgs_out_plink
mkdir -p "$OUTDIR"

declare -A BASES=(
  [ieu-b-4813]="$BASE_DIR/ieu-b-4813_base_prsice.bestP.nodup.nopal.txt.gz"
  [ieu-b-4815]="$BASE_DIR/ieu-b-4815_base_prsice.bestP.nodup.nopal.txt.gz"
  [ieu-b-4817]="$BASE_DIR/ieu-b-4817_base_prsice.bestP.nodup.nopal.txt.gz"
  [ieu-b-4845]="$BASE_DIR/ieu-b-4845_base_prsice.bestP.nodup.nopal.txt.gz"
  [ieu-b-4841]="$BASE_DIR/ieu-b-4841_base_prsice.bestP.nodup.nopal.txt.gz"
  [ieu-b-4843]="$BASE_DIR/ieu-b-4843_base_prsice.bestP.nodup.nopal.txt.gz"
  [ieu-b-4857]="$BASE_DIR/ieu-b-4857_base_prsice.bestP.nodup.nopal.txt.gz"
)

POPS=(AFR AMR EAS EUR SAS)

for trait in "${!BASES[@]}"; do
  base="${BASES[$trait]}"

  # Generate PLINK score input with 3 columns: SNP, A1, BETA
  zcat "$base" | awk 'NR>1{print $3,$4,$9}' > score.${trait}.snp_a1_beta.txt

  for pop in "${POPS[@]}"; do
    tgt="$TGT_DIR/${pop}.merged.withsex.withrs.onlyrs.nodup"
    out="$OUTDIR/${trait}.${pop}"

    # Note: In PLINK 1.9, --score column indices = SNP_col A1_col weight_col
    # Our score files have: 1=SNP, 2=A1, 3=BETA
    "$PLINK" \
      --bfile "$tgt" \
      --score score.${trait}.snp_a1_beta.txt 1 2 3 sum \
      --out "$out"

    # Output: ${out}.profile, where columns SCORESUM/AVG represent individual PGS
  done
done


# 5) Merge all .profile results by trait and population

# Step 1. Merge all .profile files into one matrix
cd /work/11063/liliu/ls6/pgs_out_plink
OUT=/work/11063/liliu/ls6/pgs_merged/PGS_matrix.tsv
mkdir -p /work/11063/liliu/ls6/pgs_merged

# Trait IDs and populations
traits=(4813 4815 4817 4841 4843 4845 4857)
pops=(AFR AMR EAS EUR SAS)

# Header
echo -e "FID\tIID\tTrait\tPop\tPGS" > "$OUT"

# Loop over all combinations
for trait in "${traits[@]}"; do
  for pop in "${pops[@]}"; do
    f="ieu-b-${trait}.${pop}.profile"
    if [[ -f "$f" ]]; then
      awk -v t="ieu-b-${trait}" -v p="$pop" 'NR>1{print $1,$2,t,p,$6}' OFS="\t" "$f" >> "$OUT"
    fi
  done
done


# Step 2. Add population (country) and super-population info
PANEL=/work/11063/liliu/ls6/gwas_20130502/integrated_call_samples_v3.20130502.ALL.panel
PGS=/work/11063/liliu/ls6/pgs_merged/PGS_matrix.tsv
OUT=/work/11063/liliu/ls6/pgs_merged/PGS_with_pop.tsv

awk 'NR==FNR{pop[$1]=$2;super[$1]=$3;next}
     NR>1{print $0, pop[$2], super[$2]}' OFS="\t" "$PANEL" "$PGS" \
> "$OUT"


# Step 3. Merge with trait names
cd /work/11063/liliu/ls6/pgs_merged

TRAIT_INFO=/work/11063/liliu/ls6/opengwas_data/opengwas_trait_id.csv
PGS=PGS_with_pop.tsv
OUT=PGS_with_pop_trait.tsv

# 1) Read the CSV file (comma-separated) as mapping: id -> trait name
awk -F',' 'NR>1{map[$1]=$2} END{for(k in map)print k"\t"map[k]}' "$TRAIT_INFO" \
  | sort -k1,1 > trait_map.tsv

# 2) Apply the mapping inside awk, output with desired column order and header
awk -v OFS="\t" 'BEGIN{
  while ((getline < "trait_map.tsv") > 0) { id=$1; name=$2; trait[id]=name }
}
NR==1{
  print "FID","IID","TraitID","TraitName","Pop","PGS","Population","SuperPop"; next
}
{
  tid=$3; tname=(tid in trait? trait[tid] : "NA")
  print $1,$2,tid,tname,$4,$5,$6,$7
}' "$PGS" > "$OUT"
