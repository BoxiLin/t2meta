# module load R/3.6.1

library(dplyr)

dat_significant_raw <- data.table::fread("dat_significant_0.01.txt")

# Table of the SNPs
urate <- dat_significant_raw %>%
  dplyr::mutate(MajorMinor = paste0( ref," / ",alt),
                maf = paste0("(", round(minor_AF,3),
                             " / ",
                            round(minor_AF.x,3),
                            " / ",
                            round(minor_AF.y,3),")")) %>%
  dplyr::select("code","variant","rsid","chr","pos","ref","alt",
              "maf", "MajorMinor", 
              "beta.x","beta.y",
              "info",
              "se.x","se.y",
              "tstat.x","tstat.y",
              "pval.x","pval.y",
              "p.T.I","p.TSG.L","p.TSG.Q") %>%
  dplyr::filter(pval.x<5e-8 & pval.y<5e-8 &
                beta.x*beta.y<0) 
urate_to_print <- urate %>%
  dplyr::select(variant,rsid, chr, pos, info, MajorMinor,# The major/minor is consistent with REF/ALT
                maf, beta.x, beta.y,
                pval.x, pval.y,
                p.T.I, p.TSG.L, p.TSG.Q) 

colnames(urate_to_print) <- c("variant","SNP", "CHR", "BP","INFO","MM","MAF",
                              "beta_{Female}", "beta_{Male}",
                              "T_{Female}", "T_{Male}", "T_{Diff}",
                              "T_{1metaL}","T_{2metaQ}")
data.table::fwrite(urate_to_print, file = "urate/s2dsnps_urate.csv")

urate_to_print <- data.table::fread(file = "urate/s2dsnps_urate.csv")