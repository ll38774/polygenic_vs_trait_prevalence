# 🧬 Polygenic vs Trait Prevalence

Assessing whether polygenic risk predictions across human populations align with observed disease and trait prevalence using within-family GWAS and 1000 Genomes data.

---

## 💡 Abstract
This project investigates whether population-level polygenic scores (PGS) correlate with the observed prevalence of complex diseases and quantitative traits across diverse human populations, using within-family GWAS and 1000 Genomes reference data.

---
## 📄 Data

1. **Base GWAS Data**  
   Within-family GWAS summary statistics from Howe *et al.*, *Nature Genetics* (2022).  
   🔗 [https://www.nature.com/articles/s41588-022-01062-7](https://www.nature.com/articles/s41588-022-01062-7)

2. **Target Data**  
   Genotype data from the **1000 Genomes Project (Phase 3, 2013 release)**.

3. **PGS Computation Tool**  
   Polygenic scores (PGS) were computed using **PLINK 1.9** with standardized SNP effect sizes (β).

4. **Trait Data Sources**

   | Trait | Source | Link |
   |:------|:--------|:-----|
   | Height | NCD-RisC (Non-Communicable Disease Risk Factor Collaboration) | 🔗 [https://ncdrisc.org/data-downloads-height.html](https://ncdrisc.org/data-downloads-height.html) |
   | Body mass index (BMI) | NCD-RisC Global BMI Database | 🔗 [https://ncdrisc.org/data-downloads-adiposity.html](https://ncdrisc.org/data-downloads-adiposity.html) |
   | Systolic blood pressure (SBP) | NCD-RisC Global Blood Pressure Database | 🔗 [https://ncdrisc.org/data-downloads-blood-pressure.html](https://ncdrisc.org/data-downloads-blood-pressure.html) |
   | LDL cholesterol | NCD-RisC Blood Lipids Database | 🔗 [https://ncdrisc.org/data-downloads-cholesterol.html](https://ncdrisc.org/data-downloads-cholesterol.html) |
   | HDL cholesterol | NCD-RisC Blood Lipids Database | 🔗 [https://ncdrisc.org/data-downloads-cholesterol.html](https://ncdrisc.org/data-downloads-cholesterol.html) |

> ⚠️ *Due to data size limitations, only trait-level summary data are provided here.*

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

---

## 📊 Results Overview

This directory contains all population-level polygenic score (PGS) results and downstream analysis outputs generated from the `polygenic_vs_trait_prevalence` project.

### Structure

- **pgs_merged/**  
  Contains the merged PGS matrix (`PGS_with_pop_trait.tsv`) combining all traits and super-populations, used for population-level correlation analyses.

- **results_bmi_pgs/**  
  Outputs of BMI-related analyses, including regression results and plots of PGS–BMI relationships across super-populations.

- **results_hdl_pgs/**  
  Results for HDL cholesterol, containing regression outputs, population scatterplots, and correlation statistics.

- **results_height_pgs/**  
  Outputs for height trait analyses, including population-level and sex-stratified correlations.

- **results_ldl_pgs/**  
  LDL cholesterol results, including cross-population comparisons and regression summaries.

- **results_sbp_pgs/**  
  Systolic blood pressure (SBP) analysis results, containing both global and population-specific association outputs.

### File Descriptions

Each subdirectory typically includes:
- `PGS_trait_summary.tsv`: Summary statistics for the PGS–trait correlation analysis.  
- `scatter_*.pdf`: Visualizations showing PGS–phenotype correlations across populations.  
- `*_regression.txt`: Regression output tables with β, r, and p-values.  

### Notes
All results were generated using:
- **Base GWAS:** Within-family GWAS (Howe et al., 2022)  
- **Target data:** 1000 Genomes Phase 3 (20130502 release)  
- **Tools:** PLINK 1.9, R (ggplot2, broom, dplyr), and custom Bash pipelines.

These outputs serve as the empirical validation for evaluating the relationship between polygenic predictions and observed population-level trait variation.

