# module load R/4.2.2

library(stringr)
library(dplyr)
library(karyoploteR)
library("latex2exp")
library(ggplot2)
library(GGally)
library(gridExtra)
library(cowplot)


###### Figure 1 (PP-Plot) ###########


maf = "_0.01"

  print(maf)
  file <- paste0("4-aggregate/dat_significant", maf, ".txt") 
  dat_significant <- data.table::fread(file) 
  
  
  my_hex <- function(data, mapping, ...){
    p <- ggplot(data = data, mapping = mapping) +
      geom_hex(bins = 50) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.5) +
      scale_fill_gradient(low = "pink", high = "darkred", 
                          limits = c(0,50), guide = FALSE) +
      coord_equal(ratio = 1) +
      expand_limits(x = range(data), y = range(data)) +
      #  theme(axis.title = element_blank())
      theme(axis.title = element_blank(),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(color = "lightgrey"),
            panel.grid.minor = element_line(color = "lightgrey"),
            panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
    p
  }
  
  
  # Define function to create diagonal labels
  my_diag <- function(data, ...){
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = colnames(data), size = 5) +
      theme_void()
  }
  
  
  # Create a colorbar; for levels only therefore can be created with arbitrary data
  colorbar <- ggplot(data=mtcars, aes(mpg, cyl)) + 
    geom_hex(bins = 50) +
    theme_void() +
    scale_fill_gradient(low = "pink", high = "darkred", 
                        limits = c(0, 50)) # adjust limits as needed
  
  # Prepare for grid plot
  plot_list <- list()
  
  
  data = as.data.frame(-log10(dat_significant[,c("pval.x","pval.y","p.T.I","p.TSG.L","p.TSG.Q")]))
  
  # Create pairwise hexbin plot
  for (i in 1:5) {
    for (j in 1:5) {
      if (j < i) {
        plot_list[[paste0(i, j)]] <- ggplotGrob(ggally_blank())
      } else if (j == i) {
        plot_list[[paste0(i, j)]] <- ggplotGrob(my_diag(data = data[,j, drop=FALSE]))
      } else {
        p <- ggally_points(data, mapping = aes_string(names(data)[j], names(data)[i]))
        plot_list[[paste0(i, j)]] <- ggplotGrob(my_hex(p$data, mapping = p$mapping))
      }
    }
  }
  
  
  # Add legend
  plot_list <- c(plot_list, list("legend" = get_legend(colorbar)))
  
  # Combine the plot matrix and the colorbar
  # Combine the plot matrix and the colorbar
  a = grid.arrange(grobs = plot_list, 
                   layout_matrix = rbind(c(1, 2, 3, 4, 5, 26),
                                         c(6, 7, 8, 9, 10, 26),
                                         c(11, 12, 13, 14, 15, 26),
                                         c(16, 17, 18, 19, 20, 26),
                                         c(21, 22, 23, 24, 25, 26)))
  
  
  cowplot::save_plot(paste0("figure1", maf,".png"), a, 
                     dpi = 800, base_width = 13, base_height = 10)
  


###### Figure 2 Manhattan ###########

snp_set <- c("T.2metaQ but not T.1metaL", "T.1metaL but not T.2metaQ",
             "T.2metaQ but not any other", "any other but not T.2metaQ")

for ( i in 1:4) {
  print(i)
if (i == 1) {
  dat_subset <- dplyr::filter(dat_significant, p.TSG.Q < 5e-8 & p.TSG.L > 5e-8)
} else if (i == 2) {
  dat_subset <- dplyr::filter(dat_significant, p.TSG.Q > 5e-8 & p.TSG.L < 5e-8)
} else if (i == 3) {
  dat_subset <- dplyr::filter(dat_significant, p.TSG.Q < 5e-8 & p.TSG.L > 5e-8 & pval.x > 5e-8 & pval.y > 5e-8 & p.T.I > 5e-8)
} else if (i == 4) {
  dat_subset <- dplyr::filter(dat_significant, p.TSG.Q > 5e-8 & (p.TSG.L < 5e-8 | pval.x < 5e-8 | pval.y < 5e-8 | p.T.I < 5e-8))
}
  

dat_granges <- makeGRangesFromDataFrame(dat_subset,start.field = "pos", end.field = "pos")
names(dat_granges) <- dat_subset$variants

png(paste("figure2-manhattan",maf,"_", snp_set[i],".png" , sep = ""),
    units = "px", res = 300,
    type = "cairo",
    height  = 5000,width = 3000)
ymax = 200
kp <- plotKaryotype(plot.type=4, labels.plotter = NULL, ideogram.plotter=NULL)
kpAddChromosomeNames(kp, col = 'black', srt = 90, cex = 1, chr.names=c(paste0("Chr ", c(1:22,"X"))," "))
kpAddMainTitle(kp, main=snp_set[i], col="red",
               cex = 1.5)


kp <- kpPlotManhattan(kp, data = dat_granges, 
                      pval = log1p(-log10(dat_subset$pval.x)),logp = FALSE, 
                      ymin = 0,ymax = 6,
                      points.cex = 0.5,
                      genomewideline = log1p(-log10(5e-8)), genomewide.col = "red",suggestiveline = 0,
                      suggestive.lty = 1,
                      r0=autotrack(5,5), points.col = "2blues")
kpAxis(kp, ymin=0, ymax=6,cex = 0.8,
       tick.pos = log1p(c(2,5,10,50,100,200)),
       labels = c(2,5,10,50,100,200),
       r0=autotrack(5,5))


kp <- kpPlotManhattan(kp, data = dat_granges,
                      pval = log1p(-log10(dat_subset$pval.y)),logp = FALSE,
                      ymin = 0,ymax = 6,
                      points.cex = 0.5,
                      genomewideline = log1p(-log10(5e-8)), genomewide.col = "red",suggestiveline = 0,
                      suggestive.lty = 1,
                      r0=autotrack(4,5), points.col = "2blues")
kpAxis(kp, ymin=0, ymax=6,cex = 0.8,
       tick.pos = log1p(c(2,5,10,50,100,200)),
       labels = c(2,5,10,50,100,200),
       r0=autotrack(4,5))

kp <- kpPlotManhattan(kp, data=dat_granges,
                      pval = log1p(-log10(dat_subset$p.T.I)),logp = FALSE,
                      points.cex = 0.5,
                      ymin = 0,ymax = 6,
                      genomewideline = log1p(-log10(5e-8)), genomewide.col = "red", suggestiveline = 0,
                      suggestive.lty = 1,
                      r0=autotrack(3,5), points.col = "2blues")
kpAxis(kp, ymin=0, ymax=6,cex = 0.8,
       tick.pos = log1p(c(2,5,10,50,100,200)),
       labels = c(2,5,10,50,100,200),
       r0=autotrack(3,5))

kp <- kpPlotManhattan(kp, data=dat_granges,
                      pval = log1p(-log10(dat_subset$p.TSG.L)),logp = FALSE,
                      points.cex = 0.5,
                      ymin=0, ymax=6,
                      genomewideline = log1p(-log10(5e-8)), genomewide.col = "red", suggestiveline = 0,
                      suggestive.lty = 1,
                      r0=autotrack(2,5), points.col = "2blues")
kpAxis(kp, ymin=0, ymax=6,cex = 0.8,
       tick.pos = log1p(c(2,5,10,50,100,200)),
       labels = c(2,5,10,50,100,200),
       r0=autotrack(2,5))

kp <- kpPlotManhattan(kp, data=dat_granges,
                      pval = log1p(-log10(dat_subset$p.TSG.Q)),logp = FALSE,
                      ymin = 0,ymax = 6,
                      points.cex = 0.5,
                      genomewideline = log1p(-log10(5e-8)), genomewide.col = "red", suggestiveline = 0,
                      suggestive.lty = 1,
                      r0=autotrack(1,5), points.col = "2blues")
kpAxis(kp, ymin=0, ymax=6,cex = 0.8,
       tick.pos = log1p(c(2, 5, 10, 50, 100, 200)),
       labels = c(2, 5, 10, 50, 100, 200),
       r0=autotrack(1,5))
dev.off()


}


