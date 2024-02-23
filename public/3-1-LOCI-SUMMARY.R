library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library("latex2exp")


mafs = c("0.05","0.01","0.01_0.05")
files = c("4-aggregate/dat_significant_0.05.txt", "4-aggregate/dat_significant_0.01.txt","4-aggregate/dat_significant_0.01_0.05.txt")

maf = 0.01
file <- paste0("dat_significant_", maf, ".txt")
signif <- data.table::fread(file)

code_list <- unique(signif$code)

#### N-N by traits plots (Loci) #####
sapply(code_list, function(code_i) {
  print(code_i)
  signif_sub <- dplyr::filter(signif, code == code_i)
  dat_x <- dplyr::filter(signif_sub, pval.x <5e-8) %>% dplyr::select(code,variant, pval.x)
  dat_y <- dplyr::filter(signif_sub, pval.y <5e-8) %>% dplyr::select(code,variant, pval.y)
  dat_i <- dplyr::filter(signif_sub, p.T.I <5e-8) %>% dplyr::select(code,variant, p.T.I)
  dat_l <- dplyr::filter(signif_sub, p.TSG.L <5e-8) %>% dplyr::select(code,variant, p.TSG.L)
  dat_q <- dplyr::filter(signif_sub, p.TSG.Q <5e-8) %>% dplyr::select(code,variant, p.TSG.Q)
  
  dd <- list(dat_x, dat_y, dat_i,dat_l,dat_q)
  index <- c("x","y","i","l","q")
  
  for (i in 1:5) {
    if (dim(dd[[i]])[1]>0) {
      out <- dd[[i]]
      names(out)[3]<- "P"
      data.table::fwrite(out, paste0("test_specific_summary/", code_i, "_", index[i],".txt"), sep = "\t")
    }
  }

})

### Run LD clumping with PLINK (script `3-2-clumping.sh`)


res <- sapply (code_list, function(code_i) {
  out <- sapply(c("x","y","i","l","q"), function(tt){
    file_path <- paste0("locus_out_100kb/",code_i,"_", tt, ".clumped")
    if(file.exists(file_path)) {
      # If the file exists, read the file and count the number of rows
      data <- data.table::fread(file_path)
      num_rows <- nrow(data)
    } else {
      # If the file does not exist, assign zero
      num_rows <- 0
    }})
  names(out) <- c("x","y","i","l","q")
  out
}) %>% t %>% as.data.frame

res$code <- rownames(res)

dict <- data.table::fread("../dict_continuous_trait_290.txt")


result <- dplyr::left_join(res, dict, by = c("code")) # %>%
#  left_join(res_new, by = "code")
result


data.table::fwrite(result, file = "result_independent_count_with_union_maf0.01_10mb.txt", sep = "\t")
