args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("wrong number", call.=FALSE)
}

code <- "30850_raw"

code <- args[1]

library(dplyr)

dat <- data.table::fread(paste("temp/",code,".txt", sep = ""))


dat$TSG.Q <- dat$tstat.y^2 + dat$tstat.x^2
dat$p.TSG.Q <- pchisq(dat$TSG.Q, df = 2, lower.tail = FALSE)

w.f = (1/dat$se.x)^2
w.m = (1/dat$se.y)^2
dat$TSG.L <- sqrt(w.f/(w.f+w.m))*dat$tstat.x+sqrt(w.m/(w.f+w.m))*dat$tstat.y
dat$p.TSG.L = 2*pnorm(abs(dat$TSG.L), lower.tail = F)

dat$T.I <- (dat$beta.x-dat$beta.y)/(sqrt(dat$se.x^2+dat$se.y^2))
dat$p.T.I <- 2*pnorm(abs(dat$T.I), lower.tail = F)

# Truncation
dat <- dplyr::mutate(dat,
                     pval.x = ifelse(dat$pval.x < 1e-200, 1e-200, dat$pval.x),
                     pval.y = ifelse(dat$pval.y < 1e-200, 1e-200, dat$pval.y),
                     pval = ifelse(dat$pval < 1e-200, 1e-200, dat$pval),
                     p.T.I = ifelse(dat$p.T.I < 1e-200, 1e-200, dat$p.T.I),
                     p.TSG.Q = ifelse(dat$p.TSG.Q < 1e-200, 1e-200, dat$p.TSG.Q),
                     p.TSG.L = ifelse(dat$p.TSG.L < 1e-200, 1e-200, dat$p.TSG.L))



data.table::fwrite(dat, file = paste("temp/",code,"-tested.txt", sep = ""), sep = "\t")


dat_significant <- dplyr::filter(dat,
                                 (pval.x < 5e-8)|
                                   (pval.y < 5e-8)|
                                   # (pval<5e-8)|
                                   (p.T.I <5e-8) |
                                   (p.TSG.Q<5e-8)|
                                   (p.TSG.L<5e-8))
data.table::fwrite(dat_significant, sep = "\t",
                   file = paste("out/significant-subset/",code,"-significant.txt",sep = ""))




