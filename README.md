# Shared genetic determinants of eGFR and albuminuria with a focus on type 2 diabetes: A large-scale European study

## Overview

This repository contains Python code used to perform the genetic association analyses described in the manuscript "Shared genetic determinants of eGFR and albuminuria with a focus on type 2 diabetes: A large-scale European study"
“Shared genetic determinants of eGFR and albuminuria with a focus on type 2 diabetes” (submitted to Kidney International).
The code implements preprocessing, regression analyses, and visualization for the investigation of shared genetic loci influencing estimated glomerular filtration rate (eGFR) and urinary albumin-creatinine ratio (UACR) in UK Biobank participants.

## Repository contents
### File	Description
article_eGFR.py	Main analysis pipeline. Loads genotype and phenotype data from UK Biobank, performs preprocessing (outlier removal, scaling, UACR and eGFR calculation), and conducts multiple linear regressions between SNPs and eGFR/UACR. Includes generation of summary statistics, scatter plots, and sensitivity analyses for both the full and type 2 diabetes (T2D) subsets.
percent_change.py	Utility script that computes and appends percentage change in mean eGFR and UACR based on regression beta estimates, producing final summary Excel files for both the overall and T2D populations.

### Data availability
Individual-level data were obtained from the UK Biobank under approved application #71699.
Due to UK Biobank data-use restrictions, raw phenotype and genotype data cannot be shared publicly.

### Requirements

Python ≥ 3.9

Packages: pandas, numpy, statsmodels, scikit-learn, matplotlib, seaborn, openpyxl, xlsxwriter, adjustText
