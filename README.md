

# Repository for Identification of suicide brain transcriptomic signatures using meta-analysis of multiple cohorts

This repository contains source code (without source data) for the analyses. Source datasets could be obtained from GEO. 
Code was written on a Linux machine and is not tested to work on Windows.

Contents:

- Data_preprocessing_analysis/TRANSCRIPT_SUICIDE_PREPR_ANALYSIS_SCRIPT.R contain primary DEs and cohort preprocessing
- Data_preprocessing_analysis/DEMOGRAPHICS_TRANSCRIPT_SCRIPT.R code for preparation of demographics
- Data_preprocessing_analysis/PrepareCountsSTAR.py, Data_preprocessing_analysis/PrepareCountsSTARControl.py, Data_preprocessing_analysis/PrepareCountsSTARprefSEP.py python files to call STAR aligner with different strategies to handle errors
- Data_preprocessing_analysis/Biomart_parser_book.ipynb Jupyter file to prepare ID mappings to gene symbols
- OLINK_PSY_SUICIDE_SCRIPT.R contains main analysis run after preprocessing
- Moderator_calculation.R contains steps to calculate cohort moderators
- Cell_expression_imputation_tests.R test runs to construct signature matrices for CIBERSORTx deconvolution
- Cell_expression_imputation_analysis.R code for cell expression deconvolution meta-analysis
- AI_gene_classification.ipynb Jupyter file to classify genes based on GPT-5