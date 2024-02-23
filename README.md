# $T_{2,metaQ}$

This document includes code to reproduce the data analysis and simulation studies in the manuscript ***"**Better together against genetic heterogeneity: a sex-combined joint main and interaction analysis of 290 quantitative traits in the UK Biobank**".***

# 0. Neale Lab’s GWAS summary statistics extraction and download

Input:

- `UKBB.csv`
    
    Downloaded from [Neale Lab’s GWAS summary statistics V3](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/) ([Documentation](https://github.com/Nealelab/UK_Biobank_GWAS))
    

Script:

- `0-PREPROCESS-PHENOTYPE.R`
    
    Extract wget commands for phenotypes for downstream analysis, stored in `0-DOWNLOAD.sh`
    
- `0-main_download_array2.sh`
    
    Download summary statistics for 290 phenotypes*3 studies (female, male and bothsex)  
    
- `0-main_download_gunzip`
    
    File conversion: *.bgz→ *.tsv
    

Output:

- `[UKB ID]_raw.gwas.imputed_v3.{female, male, both_sexes}.tsv`
    
    Summary statistics across 290 phenotypes, each with three studies (female, male, both_sexes)
    
- `dict_continuous_trait_290.txt`
    
    A two-column (UKB ID, Trait name) dictionary for phenotype analyzed.
    

# 1. QC, run tests and store results

Input:

- `[UKB ID]_raw.gwas.imputed_v3.{female, male, both_sexes}.tsv`
- `variants_updated.tsv`
    
    Variant annotation files downloaded from  [Neale Lab’s GWAS summary statistics V3](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/) ([Documentation](https://github.com/Nealelab/UK_Biobank_GWAS))
    

Script:

- `1-main_analysis.sh`
    - `1-1-QC.R`
        
        Running QC, and extract relevant columns
        
    - `1-2-TEST.R`
        
        Calculate the five association tests discussed in manuscript, extract signals
        
    - `1-4-AGGREGATE-SUMMARY.R`
        
        Aggregate signals across all 290 phenotypes; SNP annotation
        

Output:

- `dat_significant_0.01.txt`
    
    SNP-Phenotype associated with genome-wide significance in at least one of the five tests considered
    

# 2. Plot Figures 1 (pairwise-scatter) and 2 (stacked-manhattan)

Input:

- `dat_significant_0.01.txt`

Script:

- `2-PLOT-FIGURE1-2.R`
    
    Create Figures 1 and 2.
    

Output

- `figure1_maf_0.01.png`
    
    Figure 1.
    
- `figure2-manhattan_maf_0.01_{T.2metaQ but not T.1metaL, T.1metaL but not T.2metaQ}.png`
    
    Sub-panels of Figure 2. Testosterone-removed results are available by filtering out testosterone results (UKB code: 30850).
    

# 3. Plot Figure 3 (by-trait comparisons for SNPs and loci counts respectively)

Input:

- `dat_significant_0.01.txt`
- `dict_continuous_trait_290.txt`
- `EUR_phase3_nodup.{bed,bim,fam}`
    
    QCed UK Biobank genotype data as reference for LD clumping
    

Script:

- `3-1-LOCI-SUMMARY.R` (1)
    
    Extract trait-and-test-specific summary statistics for LD clumping
    
    - `3-2-clumping.sh`
        
        LD clumping with PLINK 1.90
        
- `3-1-LOCI-SUMMARY.R` (2)
    
    Summarize LD clumping results.
    
- `3-4-PLOT-FIGURES.R`
    
    Create Figure 3.
    

Output:

- `figure 3-snps_loci_10mb.png`
    
    Figure 3.
    

# 4.  Create Table 2 (testosterone signals)

Input:

- `dat_significant_0.01.txt`

Script:

- `5-FUMA-PROCESS.R` (1)
    
    Extract SNPs uniquely identified from $T_{2,metaQ}$
    
    - `FUMA GWAS SNP2GENE`
        
        Use web-based tool [FUMA](https://fuma.ctglab.nl/) (Watanabe et al., 2017) to characterize genomic loci and annotate candidate SNPs in genomic loci based on GWAS Catalog ([Documentation](https://fuma.ctglab.nl/tutorial#snp2gene)).
        
- `5-FUMA-PROCESS.R` (2)
    - Process the FUMA output files.

Output:

- `loci_36_testosterone.csv`
    
    Raw data to create Table 2.
    

# 5.  Create Table 3 (urate signals)

Input:

- `dat_significant_0.01.txt`
    
    Annotation for testosterone loci uniquely identified from $T_{2,metaQ}$
    

Script:

- `5-URATE-SNPS-W-OPPOSITE-EFFECT.R`
    
    Urate signals with genome-wide significance in both sexes but in opposite directions
    

Output:

- `s2dsnps_urate.csv`
    
    Raw data to create Table 3.
    

# Appendix S2, S3: Simulation studies

`cd simulation/`

Script:

- `simulation_library.R`
    
    Helper function for running the five association tests; parameter grid generators
    
- `simulation_t1e.R`
    
    Run type I error evaluation
    
- `simulation_power.R`
    
    Run power evaluation at genome-wide significant level
    
- `combine.R`

Output:

- Figures in Appendix S2 and Appendix S3
