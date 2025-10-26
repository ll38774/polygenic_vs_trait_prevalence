#!/usr/bin/env bash
# ===============================================
# Script Name: 01_convert_vcf_to_prsice.sh
# Description: Convert GWAS-VCF files to PRSice-compatible base files.
# Author: Li Liu
# Usage: bash 01_convert_vcf_to_prsice.sh
# ===============================================

# Preparation base data:
# Download VCF data from OpenGWAS manually before running this script.

# Convert VCF Files to PRSice-Compatible Base Files
for f in *.vcf.gz; do
  prefix=${f%.vcf.gz}
  echo "Processing $f ..."

  (
    echo -e "CHR\tBP\tSNP\tA1\tA2\tN\tSE\tP\tBETA\tINFO\tMAF"
    zcat "$f" | awk -F'\t' -v OFS='\t' '
      /^##/ {next}
      /^#/  {next}
      {
        chr=$1; bp=$2; snp=$3; ref=$4; alt=$5;
        split($9, fk,  ":");
        split($10, fv, ":");

        n=""; se=""; lp=""; es=""; si=""; af="";
        for(i=1;i<=length(fk);i++){
          key=fk[i]; val=fv[i];
          if(key=="SS"){n=val}
          else if(key=="SE"){se=val}
          else if(key=="LP"){lp=val}
          else if(key=="ES"){es=val}
          else if(key=="SI"){si=val}
          else if(key=="AF"){af=val}
        }
        p = (lp=="" ? "" : exp(-lp*log(10)))
        if(es!=""){ print chr, bp, snp, alt, ref, n, se, p, es, si, af }
      }'
  ) > "${prefix}_base_prsice.txt"

  echo "Finished: ${prefix}_base_prsice.txt"
done

# Verify Correct Generation of PRSice Base Files
# Display the first 5 lines of each generated file in a formatted table
head -5 *_base_prsice.txt | column -t

# Check the number of columns (should print 11)
awk 'NR>1{print NF; exit}' *_base_prsice.txt   # should print 11
