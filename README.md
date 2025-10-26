# 🧬 Polygenic vs Trait Prevalence

Assessing whether polygenic risk predictions across human populations align with observed disease and trait prevalence using within-family GWAS and 1000 Genomes data.

---

## 📄 Abstract
This project investigates whether population-level polygenic scores (PGS) correlate with the observed prevalence of complex diseases and quantitative traits across diverse human populations, using within-family GWAS and 1000 Genomes reference data.

---

## ⚙️ Code Overview
This repository contains a complete pipeline for computing and analyzing population-level PGS from GWAS summary statistics and 1000 Genomes genotype data.  
All scripts are written for **PLINK 1.9** and **LDSC**, designed to run on HPC systems (e.g., TACC).

| No. | Script |
|-----|---------|
| **01_convert_vcf_to_prsice.sh** | Convert GWAS-VCF files from OpenGWAS into PRSice-compatible base files using PLINK 1.9. |
| **02_ldsc_heritability_qc.sh** | Perform heritability QC (*h²ₛₙₚ > 0.05*) for base GWAS data using LDSC on PLINK 1.9-formatted files. |
| **03_deduplicate_snps.sh** | Remove duplicate SNP entries and keep variants with the most significant (lowest) *P*-value. |
| **04_remove_ambiguous_snps.sh** | Filter out non-ACGT and palindromic (A/T, T/A, C/G, G/C) SNPs to ensure allele-strand consistency. |
| **05_build_plink_per_superpop.sh** | Build population-specific PLINK datasets from 1000 Genomes VCFs (EUR/AFR/EAS/SAS/AMR) using PLINK 1.9. |
| **06_map_rsid_and_dedup.sh** | Map rsIDs from base GWAS to target PLINK datasets, verify overlap, retain rsID-only SNPs, and remove duplicates. |
| **07_pgs_calculation_and_merge.sh** | Calculate PGS for multiple traits and populations using PLINK 1.9, then merge all outputs with population and trait annotations. |
| **\*_correlation.R** | Perform correlation analyses between PGS and population-level traits for each quantitative phenotype. |
