#!/usr/bin/env bash
# ===============================================================
# Script Name: 06_map_rsid_and_dedup.sh
# Description: Map rsIDs from base GWAS data to target .bim files,
#              verify overlap coverage, align PLINK file triplets,
#              keep only rsID SNPs, and remove duplicate rsIDs.
# NOTE: Run separately for each super-population (EUR/AFR/EAS/SAS/AMR).
# ===============================================================


# rsID mapping
cd /work/11063/liliu/ls6/gwas_20130502
PLINK=/work/11063/liliu/ls6/tools/plink1.9/plink
BASE=/work/11063/liliu/ls6/opengwas_data/sibling_gwas/base_gwas_h2_0.05/ieu-b-4813_base_prsice.txt

# 1) Extract mapping table from base: CHR BP A1 A2 rsID
#    If CHR column in base file contains "chr" prefix, use:
#    awk 'NR>1{gsub(/^chr/,"",$1); print $1,$2,$4,$5,$3}'
awk 'NR>1{print $1,$2,$4,$5,$3}' "$BASE" > base_chr_bp_a1_a2_rsid.txt

# 2) Replace column 2 in target .bim with rsID (consider A1/A2 swap + complement chain)
BIM_IN=qc_AMR/plink/AMR.merged.withsex.bim
BIM_OUT=qc_AMR/plink/AMR.merged.withsex.withrs.bim

awk 'BEGIN{OFS="\t"}
function comp(b){return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":b}
FNR==NR{
  # Build multiple key mappings from base -> rsID
  k1=$1":"$2":"$3":"$4;                 map[k1]=$5
  k2=$1":"$2":"$4":"$3;      if(!(k2 in map))  map[k2]=$5
  kc1=$1":"$2":"comp($3)":"comp($4);  if(!(kc1 in map)) map[kc1]=$5
  kc2=$1":"$2":"comp($4)":"comp($3);  if(!(kc2 in map)) map[kc2]=$5
  next
}
{
  key=$1":"$4":"$5":"$6
  if(key in map){ $2=map[key] }   # Replace SNP column with rsID when matched
  print
}' base_chr_bp_a1_a2_rsid.txt "$BIM_IN" > "$BIM_OUT"

# 3) Count rsID mapping success rate
awk '{if($2 ~ /^rs/) rs++; else other++} END{print "AMR rsID:",rs," non-rs:",other," total:",NR}' "$BIM_OUT"


# Check overlap and coverage with base data
# Adjust paths accordingly
BASE=/work/11063/liliu/ls6/opengwas_data/sibling_gwas/base_gwas_h2_0.05/ieu-b-4813_base_prsice.txt
BIM=/work/11063/liliu/ls6/gwas_20130502/qc_AFR/plink/AFR.merged.withsex.withrs.bim

# 1) Unique rsID count in base (~7.17 million expected)
awk 'NR>1{print $3}' "$BASE" | sort -u > base.rsids.txt
wc -l base.rsids.txt

# 2) Unique rsID count in target data (after mapping)
awk '$2~/^rs/{print $2}' "$BIM" | sort -u > target.rsids.txt
wc -l target.rsids.txt

# 3) Intersection of both (maximum number of scoreable SNPs)
comm -12 base.rsids.txt target.rsids.txt | wc -l

# 4) Calculate coverage (intersection / base)
OVL=$(comm -12 base.rsids.txt target.rsids.txt | wc -l)
TOT=$(wc -l < base.rsids.txt)
awk -v a=$OVL -v b=$TOT 'BEGIN{printf("Coverage vs BASE: %d / %d = %.2f%%\n", a,b,100*a/b)}'

# Example expected output:
# Coverage vs BASE: 7145063 / 7170322 = 99.65%


# 1) Align PLINK triplets (use .bim with rsIDs)
#    Note: Change population name accordingly (run inside /work/11063/liliu/ls6/gwas_20130502/qc_(POP)/plink)
#    Generate a triplet with prefix SAS.merged.withsex.withrs
cp SAS.merged.withsex.bed SAS.merged.withsex.withrs.bed
cp SAS.merged.withsex.fam SAS.merged.withsex.withrs.fam

# Quick check that the PLINK triplet is complete
ls -lh SAS.merged.withsex.withrs.{bed,bim,fam}


# 2) (Recommended) Keep only SNPs with rsIDs and generate an “onlyrs” version
# Extract rsID list
awk '$2 ~ /^rs/{print $2}' SAS.merged.withsex.withrs.bim > SAS.rs.keep.txt

# Filter rsID SNPs and write a clean prefix
/work/11063/liliu/ls6/tools/plink1.9/plink \
  --bfile SAS.merged.withsex.withrs \
  --extract SAS.rs.keep.txt \
  --make-bed \
  --out SAS.merged.withsex.withrs.onlyrs

# Verify output
ls -lh SAS.merged.withsex.withrs.onlyrs.{bed,bim,fam}

# Check SNP count after filtering (only rsIDs)
wc -l SAS.merged.withsex.withrs.onlyrs.bim


# Remove duplicated rsIDs
cd /work/11063/liliu/ls6/gwas_20130502/qc/plink

# 1) Use awk to generate a BIM with unique IDs and record duplicates (2nd, 3rd, etc.)
awk -v OFS="\t" '
{
  id=$2
  cnt[id]++
  if (cnt[id]==1) {
    # 1st occurrence: keep original ID
    print $1,$2,$3,$4,$5,$6 >> "SAS.tmp.uniq.bim"
  } else {
    # 2nd+ occurrence: add suffix __dupN and record in exclusion list
    new_id = id "__dup" cnt[id]
    print $1,new_id,$3,$4,$5,$6 >> "SAS.tmp.uniq.bim"
    print new_id >> "exclude_dups.txt"
  }
}' SAS.merged.withsex.withrs.onlyrs.bim

# 2) Use PLINK to read “original bed/fam + new bim” and create unified dataset
/work/11063/liliu/ls6/tools/plink1.9/plink \
  --bed SAS.merged.withsex.withrs.onlyrs.bed \
  --bim SAS.tmp.uniq.bim \
  --fam SAS.merged.withsex.withrs.onlyrs.fam \
  --make-bed \
  --out SAS.tmp.uniq

# 3) Exclude all __dupN entries (duplicates beyond the first) to get final non-duplicated dataset
/work/11063/liliu/ls6/tools/plink1.9/plink \
  --bfile SAS.tmp.uniq \
  --exclude exclude_dups.txt \
  --make-bed \
  --out SAS.merged.withsex.withrs.onlyrs.nodup

# Should be 0 (no duplicate rsIDs remain)
awk '{print $2}' SAS.merged.withsex.withrs.onlyrs.nodup.bim | sort | uniq -d | wc -l
wc -l SAS.merged.withsex.withrs.onlyrs.nodup.bim

# 4) Validate with PLINK (if bed/bim/fam mismatch, it will error)
 /work/11063/liliu/ls6/tools/plink1.9/plink \
  --bfile SAS.merged.withsex.withrs.onlyrs.nodup \
  --freq \
  --out nodup_check
