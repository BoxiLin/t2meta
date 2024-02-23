library("dplyr")
library(stringr)


dat_significant_30850 <- data.table::fread("dat_significant_0.01.txt")  %>%
  dplyr::filter(code == "30850_raw" & p.TSG.L>5e-8 & p.TSG.Q < 5e-8 & 
                  pval.x  >5e-8  & pval.y >5e-8 & p.T.I >5e-8)

fwrite(dat_significant_30850, "nealelab-testosterone-t2metaq-unique.txt", sep = "\t")


### Submit to FUMA (https://fuma.ctglab.nl/) with "nealelab-testosterone-t2metaq-unique.txt"

loci <- data.table::fread("FUMA_job252346/GenomicRiskLoci.txt")

gwas_summary <- data.table::fread("nealelab-testosterone-t2metaq-unique.txt") %>% 
  dplyr::mutate(chr = str_replace_all(chr, "chr", ""))
gwas_summary$chr[gwas_summary$chr == "X"] <- "23"
gwas_summary$chr <- as.numeric(gwas_summary$chr)


gwas_catalog <- data.table::fread("FUMA_job252346/gwascatalog.txt") 

loci_catalog = left_join(
          gwas_catalog[,c("IndSigSNP","chr","bp","GenomicLocus","snp", 
                                                      "Date", "FirstAuth","Trait", "InitialN","MappedGene","PMID",
                                                      "Context")], 
          loci[,c("GenomicLocus", "rsID",  "uniqID")],
          by = "GenomicLocus") %>%
  dplyr::filter(rsID==snp)


data.table::fwrite(loci_catalog, file = "loci_full_testosterone.csv")

dat <- data.table::fread(file = "loci_full_testosterone.csv")
dat <- unique(dat)


data.table::fwrite(dat, file = "loci_full_testosterone.csv")

loci_catalog <- loci_catalog %>%
  mutate(trait.testosterone = ifelse( rsID==snp,ifelse(str_detect(Trait, "stostero|hormone"), Trait, NA), NA))


loci_summary <- left_join(loci[,c("GenomicLocus", "rsID",  "uniqID","chr","pos")], 
                          gwas_summary[,c("variant","chr","bp", "beta.x","se.x","beta.y", "se.y","pval.x", "pval.y", "p.T.I","p.TSG.L","p.TSG.Q")],
          by = c( "chr" ="chr", "pos"="bp")) %>%
  mutate(minor_allele = sub("^[^:]*:[^:]*:[^:]*:", "", variant))

B1 <- sapply(stringr::str_extract_all(loci_summary$uniqID, "(?<=:)[^:]+(?=:)"), `[`, 2)
B2 <- stringr::str_extract(loci_summary$uniqID, "(?<=:)[^:]+$")
loci_summary$major_allele <-ifelse(loci_summary$minor_allele==B1, B2, ifelse(loci_summary$minor_allele==B2, B1, "CHECK"))
data.table::fwrite(loci_summary, file = "loci_36_testosterone.csv")