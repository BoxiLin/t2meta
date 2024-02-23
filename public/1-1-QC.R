args = commandArgs(trailingOnly=TRUE)

code <- "30850_raw"
code = args[1]

female_raw <- paste0(code, ".gwas.imputed_v3.female.tsv")
male_raw <-  paste0(code, ".gwas.imputed_v3.male.tsv")
both_raw <-  paste0(code, ".gwas.imputed_v3.both_sexes.tsv")

dat_annote <- data.table::fread("neale_lab_documents/variants_updated.tsv")

if (length(args)!=1) {
  stop("wrong number", call.=FALSE)
}

library(dplyr)

# 1. QC

dat_female_raw <- data.table::fread(paste("temp/",female_raw, sep = ""))
dat_male_raw <- data.table::fread(paste("temp/",male_raw, sep = ""))
dat_bothsex_raw <- data.table::fread(paste("temp/",both_raw, sep = ""))

dat_raw <- dplyr::left_join(dat_annote, dat_female_raw,  by = "variant") %>%
  dplyr::inner_join(dat_male_raw,  by = "variant") %>%
  dplyr::inner_join(dat_bothsex_raw, by = "variant")


dat_raw_updated <- dat_raw 
# %>% 
#   dplyr::mutate(beta.x = ifelse(minor_allele.x != minor_allele_meta, -beta.x, beta.x),
#                 tstat.x = ifelse(minor_allele.x != minor_allele_meta, -tstat.x, tstat.x),
#                 beta.y = ifelse(minor_allele.y != minor_allele_meta, -beta.y, beta.y),
#                 tstat.y = ifelse(minor_allele.y != minor_allele_meta, -tstat.y, tstat.y),
#                 beta = ifelse(minor_allele != minor_allele_meta, -beta, beta),
#                 tstat = ifelse(minor_allele != minor_allele_meta, -tstat, tstat))




dat <- dplyr::filter(dat_raw_updated, 
                    # (minor_AF.x>0.01) & (minor_AF.y>0.01) &
                    #  (low_confidence_variant.x==FALSE) & (low_confidence_variant.y==FALSE) &
                       (!is.na(beta.x)) & (!is.na(beta.y)) & (!is.na(beta)) &
                       (!is.na(se.x)) & (!is.na(se.y)) & (!is.na(se)) &
                       (!is.na(tstat.x)) & (!is.na(tstat.y)) & (!is.na(tstat)) &
                       (!is.na(pval.x)) & (!is.na(pval.y)) & (!is.na(pval))) %>%
  dplyr::select(-c(# QC criteria from meta-file
                    "info", "call_rate", "AC_meta",
                   "AF", "minor_AF_meta",     "p_hwe", "n_called", "n_not_called",
                   "n_hom_ref", "n_het", "n_hom_var", "n_non_ref", "r_heterozygosity", "r_het_hom_var",
                    "r_expected_het_frequency",
                   
                   # Annotation from meta-file
                    "ref",  "alt"  , "varid", "consequence" ,
                   "consequence_category", 
                   
                   # Available from Phenotype_summary
                   "n_complete_samples.x", "n_complete_samples.y","n_complete_samples",
                   
                   # Not aplicable after sign fliping
                   "ytx.x","ytx.y",   "ytx"))

data.table::fwrite(dat, file = paste("temp/", code,".txt",sep = ""), sep = "\t")

