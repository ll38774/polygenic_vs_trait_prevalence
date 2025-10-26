#!/usr/bin/env bash
# ===============================================================
# Script Name: build_plink_per_superpop.sh
# Description: Build per-super-population PLINK dataset from 1000G VCFs,
#              clean and deduplicate IDs (special handling for chr 8/12/14/17),
#              merge chromosomes, and add sex information.
# NOTE: Run this script separately for each super-population (EUR, AFR, EAS, SAS, AMR)
# ===============================================================

# NOTE: This pipeline is executed separately for each of the five super-populations.
# Remember to replace the population name (e.g., EUR) accordingly in the code below.

# Working directory
cd /work/11063/liliu/ls6/gwas_20130502

# Path to PLINK executable
PLINK=/work/11063/liliu/ls6/tools/plink1.9/plink

# Output directories
mkdir -p qc/tmp qc/plink

# Step 0: Extract target population (EUR) sample list
# In the 2013 1000G panel, the third column represents the super-population (EUR/AFR/EAS/SAS/AMR)
awk '$3=="EUR"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > EUR.samples.txt

# Step 1: Convert VCF to PLINK binaries (per chromosome)
# For each chromosome, extract A/C/G/T SNPs for the EUR population and generate .EUR.{bed,bim,fam}
for CHR in {1..22}; do
  VCF=ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
  OUT=qc/tmp/chr${CHR}.EUR
  $PLINK \
    --vcf $VCF \
    --double-id \
    --keep EUR.samples.txt \
    --snps-only just-acgt \
    --make-bed \
    --out $OUT
done

# Step 2: Clean, deduplicate, and assign IDs to generate final .EUR.bi files
# A. Normal chromosomes (excluding 8, 12, 14, and 17)
# These chromosomes can use --set-missing-var-ids to generate unique IDs after a single deduplication step.
for CHR in 1 2 3 4 5 6 7 9 10 11 13 15 16 18 19 20 21 22; do
  IN=qc/tmp/chr${CHR}.EUR

  # 1) Basic cleaning (remove non-biallelic SNPs; keep original IDs)
  $PLINK --bfile ${IN} --biallelic-only strict --make-bed --out ${IN}.tmp

  # 2) Keep the first occurrence per key “chr:pos:a1:a2” (TAB-delimited range file)
  awk 'BEGIN{OFS="\t"}{
    k=$1":"$4":"$5":"$6; if(!(k in s)){print $1,$4,$4,"keep"NR; s[k]=1}
  }' ${IN}.tmp.bim > ${IN}.tmp.keep.range

  $PLINK --bfile ${IN}.tmp --extract range ${IN}.tmp.keep.range \
         --make-bed --out ${IN}.tmp2

  # 3) Use PLINK to assign unique IDs (works for normal chromosomes)
  $PLINK --bfile ${IN}.tmp2 \
         --set-missing-var-ids @:#:\$1:\$2 \
         --make-bed \
         --out qc/tmp/chr${CHR}.EUR.bi
done

# B. Special chromosomes 8, 12, 14, and 17 (require special handling)
# These chromosomes contain many ID='.' entries and repeated chr:pos:A1:A2 keys;
# --set-missing-var-ids will fail due to duplicate positions.
# Apply two rounds of “keep-first” filtering (tmp → tmp2 → tmp3)
# and then directly modify the .bim ID column (column 2) to a unique string (chr:pos:a1:a2_lineNumber).
for CHR in 8 12 14 17; do
  IN=qc/tmp/chr${CHR}.EUR

  # 1) Basic cleaning
  $PLINK --bfile ${IN} --biallelic-only strict --make-bed --out ${IN}.tmp

  # 2) First “keep-first” round (tmp -> tmp2)
  awk 'BEGIN{OFS="\t"}{
    k=$1":"$4":"$5":"$6; if(!(k in s)){print $1,$4,$4,"k"NR; s[k]=1}
  }' ${IN}.tmp.bim > ${IN}.tmp.keep.range

  $PLINK --bfile ${IN}.tmp --extract range ${IN}.tmp.keep.range \
         --make-bed --out ${IN}.tmp2

  # 3) Second “keep-first” round (tmp2 -> tmp3)
  awk 'BEGIN{OFS="\t"}{
    k=$1":"$4":"$5":"$6; if(!(k in s)){print $1,$4,$4,"k"NR; s[k]=1}
  }' ${IN}.tmp2.bim > ${IN}.tmp2.keep.range

  $PLINK --bfile ${IN}.tmp2 --extract range ${IN}.tmp2.keep.range \
         --make-bed --out ${IN}.tmp3

  # 4) Modify column 2 (variant ID) in .bim to a unique ID format “chr:pos:a1:a2_lineNumber”
  cp ${IN}.tmp3.bim ${IN}.tmp3.bim.bak
  awk 'BEGIN{OFS="\t"}{
    id=$1":"$4":"$5":"$6"_"NR; $2=id; print
  }' ${IN}.tmp3.bim > ${IN}.tmp3.bim.uniq
  mv ${IN}.tmp3.bim.uniq ${IN}.tmp3.bim

  # 5) Copy as final output .EUR.bi
  cp ${IN}.tmp3.bed qc/tmp/chr${CHR}.EUR.bi.bed
  cp ${IN}.tmp3.bim qc/tmp/chr${CHR}.EUR.bi.bim
  cp ${IN}.tmp3.fam qc/tmp/chr${CHR}.EUR.bi.fam
done


# Merge step: must be executed on a compute node
cd /work/11063/liliu/ls6/gwas_20130502
PLINK=/work/11063/liliu/ls6/tools/plink1.9/plink

# Generate merge list (chr2..22)
> qc/tmp/merge_list.bi.txt
for CHR in {2..22}; do echo qc/tmp/chr${CHR}.EUR.bi >> qc/tmp/merge_list.bi.txt; done

# Perform merge (with sufficient memory)
$PLINK \
  --bfile qc/tmp/chr1.EUR.bi \
  --merge-list qc/tmp/merge_list.bi.txt \
  --make-bed \
  --memory 200000 \
  --out qc/plink/EUR.merged


# Add sex information
cd /work/11063/liliu/ls6/gwas_20130502
PLINK=/work/11063/liliu/ls6/tools/plink1.9/plink

# 1) Generate EUR sex mapping file: FID IID SEX (SEX: 1=male, 2=female, 0=missing)
awk 'NR>1 && $3=="EUR"{
  g=tolower($4);
  if     (g ~ /^m(ale)?$/)    sex=1;
  else if(g ~ /^f(emale)?$/)  sex=2;
  else                        sex=0;
  print $1, $1, sex
}' integrated_call_samples_v3.20130502.ALL.panel > EUR.sex.txt

# 2) Update the merged dataset (FID=IID=sample name, so they match directly)
$PLINK \
  --bfile qc/plink/EUR.merged \
  --update-sex EUR.sex.txt \
  --make-bed \
  --out qc/plink/EUR.merged.withsex

# 3) Verify sex counts in .fam
awk '{c[$6]++} END{print "fam male(1):",c[1]; print "fam female(2):",c[2]; print "fam missing(0):",c[0]; print "total:",NR}' \
  qc/plink/EUR.merged.withsex.fam
