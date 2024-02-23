# Aggregate significant files;

# module load R/3.6.1
# module load zlib


require(Cairo)
library(stringr)
library(dplyr)
library(karyoploteR)
library(latex2exp)
library(ggplot2)

f_list <- list.files("out/significant-subset/",full.names = TRUE) 
dat_significant <- f_list %>% 
  lapply(function(f) {
    print(f)
    code <- stringr::str_extract(f, "(?<=//)(.*)(?=-significant)")
    df.sig <-  data.table::fread(f) %>% dplyr::mutate(chr = as.character(chr))
    n = dim(df.sig)[1]
    if (n>0) {
      return(data.frame(code = code, df.sig))
    } else return(NULL)
  }) %>% 
  bind_rows 


snp_annotate <- data.table::fread("neale_lab_documents/variants_updated.tsv")

dat_significant_annotated <- left_join(dat_significant, snp_annotate[,c("variant", "info", "call_rate", "AC_meta",
                           "AF", "minor_AF_meta",     "p_hwe", "n_called", "n_not_called",
                           "n_hom_ref", "n_het", "n_hom_var", "n_non_ref", "r_heterozygosity", "r_het_hom_var",
                           "r_expected_het_frequency",
                           "ref",  "alt"  , "varid", "consequence" ,
                           "consequence_category")],
          by = "variant")

data.table::fwrite(dat_significant_annotated, file = "4-aggregate/significant_signals.txt", 
                   sep = "\t")


dat_significant_annotated <- data.table::fread("4-aggregate/significant_signals.txt")


dat_significant_0.001_0.01 <- dplyr::filter(dat_significant_annotated, !(minor_AF.x > 0.01 & minor_AF.y > 0.01))
dat_significant_0.01 <- dplyr::filter(dat_significant_annotated, minor_AF.x > 0.01 & minor_AF.y > 0.01)
dat_significant_0.05 <- dplyr::filter(dat_significant_annotated, minor_AF.x > 0.05 & minor_AF.y > 0.05)
dat_significant_0.01_0.05 <- dplyr::filter(dat_significant_0.01, !((minor_AF.x > 0.05) & (minor_AF.y > 0.05)))

data.table::fwrite(dat_significant_0.01, file = "4-aggregate/dat_significant_0.01.txt", 
       sep = "\t")


data.table::fwrite(dat_significant_0.05, file = "4-aggregate/dat_significant_0.05.txt", 
       sep = "\t")

data.table::fwrite(dat_significant_0.01_0.05, file = "4-aggregate/dat_significant_0.01_0.05.txt", 
                   sep = "\t")


data.table::fwrite(dat_significant_0.001_0.01, file = "4-aggregate/dat_significant_0.001_0.01.txt", 
                   sep = "\t")

