library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library("latex2exp")

maf <- "_0.01"
dict <- data.table::fread("../dict_continuous_trait_290.txt")

# Make plot start here

print(maf)
file <- paste0("dat_significant", maf, ".txt") 
signif <- data.table::fread(file)

code_list <- unique(signif$code)

result <- lapply(code_list, function(code_i) {
  print(code_i)
  signif_sub <- dplyr::filter(signif, code == code_i)
  count = colSums(signif_sub[,c("pval.x", "pval.y", "p.T.I", "p.TSG.L", "p.TSG.Q")]<5e-8)
  n = dim(signif_sub)[1]
  lq = sum(signif_sub$p.TSG.L < 5E-8 | signif_sub$p.TSG.Q <5E-8)
  data.frame("code" = code_i, data.frame(t(count)), n, lq)
}) %>% bind_rows() %>% left_join(dict, by = "code") 

result_loci  <- data.table::fread("result_independent_count_with_union_maf0.01_10mb.txt")



library(scales)

fancy_scientific <- function(l) {
  sapply(l, function(x) {
    if(is.na(x) || is.infinite(x) || x == 0) return(0) 
    else return(as.expression(bquote(10^.(log(x, 10)))))
  })
}

format_exp <- function(x) {
  sapply(x, function(y) {
    # If the value is NA, return NA
    if (is.na(y)) return(NA)
    
    # If the value is 0, return the expression for 0
    if (y == 0) return(expression(0))
    
    # Convert the number to scientific notation
    sci_y <- formatC(y, format = "e", digits = 1)
    
    # Split the scientific notation at the "e"
    parts <- unlist(strsplit(sci_y, "e"))
    
    mantissa <- as.numeric(parts[1])
    exponent <- as.numeric(parts[2])
    
    # Return the expression with the desired formatting
    parse(text = paste(mantissa, " %*% 10^{", exponent, "}", sep = ""))
  })
}

p1 =  ggplot(result, aes(p.TSG.Q,p.TSG.L, label = code)) +
  labs(y = TeX("Traditional meta-analysis (${T_{1, metaL}}$)"), x = TeX("Omnibus meta-analysis (${T_{2, metaQ}}$)"))+
  geom_point(colour = "darkgrey")+
  geom_text_repel(data  =result, 
                  aes(label = trait), 
                  size = 4)+
  geom_text_repel(data =  subset(result, code == "30850_raw"), aes(label = trait), size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed")  +
  scale_x_continuous() +
  scale_y_continuous()

p2 =  ggplot(result_loci, aes(q,l, label = code)) +
  labs(y = TeX("Traditional meta-analysis (${T_{1, metaL}}$)"), x = TeX("Omnibus meta-analysis (${T_{2, metaQ}}$)"))+
  geom_point(colour = "darkgrey")+
  geom_text_repel(data  = result_loci, aes(label = trait), size = 4) +
  geom_text_repel(data =  subset(result_loci, code %in% c("30850_raw","30290_raw","21021_raw")), aes(label = trait), size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed")  +
  scale_x_continuous() +
  scale_y_continuous()

pp <- plot_grid(p1, p2, ncol=2,
                labels = c("(a)", "(b)"), label_size = 16, 
                label_fontfamily = "sans", 
              #  label_fontface = "bold",
                hjust = -0.1, vjust = 1.1)
ggsave(paste0("figure3-snps_loci_10mb.png"), plot = pp, width = 14, height =7 , dpi = 300)

pp



