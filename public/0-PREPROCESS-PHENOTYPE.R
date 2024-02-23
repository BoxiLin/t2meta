### Extract phenotype for downstream analysis ##########
### Prepare meta-file summarizing each study ##########


### Extract phenotype for downstream analysis ##########
dat_raw <- data.table::fread("UKBB.csv")


dat <- dat_raw[grep("raw", dat_raw$`Phenotype Code`),c("Phenotype Code",
                                                       "Phenotype Description",
                                                       "Sex",
                                                       "wget command")]

colnames(dat) <- c("code","trait","sex","wget")

# Exclude studies with sex!=3
exclude_phenotypes <- names(which(table(dat$code)!=3))
dat_sex3 <- dat[!(dat$code %in% exclude_phenotypes),]
dat_sex3$id <- cumsum(!duplicated(dat_sex3$code))


# Unified the outputfiles
for (i in 1:dim(dat_sex3)[1]) {
  if (grepl("varorder",dat_sex3$wget[i] , fixed = TRUE)) {
    raw = stringr::str_split(dat_sex3$wget[i], " ")[[1]]
    dat_sex3$wget[i] <- paste(raw[1],raw[2],raw[3], gsub("varorder.","",raw[4]), sep = " ")
  }
}

data.table::fwrite(dat_sex3, file = "temp/0-neale-list-raw.txt", sep = "\t")


dat_sex3 = data.table::fread("temp/0-neale-list-raw.txt")

write.table(dat_sex3$wget, file = "0-DOWNLOAD.sh", col.names = F, row.names = F,quote = F)

# # Check the list of traits missed due to missed sex analysis.:
# library(dplyr)
# library(stringr)
# dat_missed <- dat[(dat$code %in% exclude_phenotypes),1:3]
# dat_missed_merged <- dat_missed %>%
#   group_by(code, trait) %>%
#   summarise(sex = str_c(sex, collapse = ", ")) %>% ungroup
# data.table::fwrite(dat_missed_merged, "out/0-phenotype/phenotype_removed_due_to_missed_sex_analysis.txt", sep = "\t")
