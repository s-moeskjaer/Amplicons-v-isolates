setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/divergence/')
library(ggplot2)
library(dplyr)
library(multcompView)
library(agricolae)


envir <- read.table('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/divergence/zz_sample_overview.csv', header = T, sep = ";")

#==================================================================================================================#
#                                                                                                                  #
#                                           NODD2+NODA2 FILTERED                                                   #
#                                                                                                                  #
#==================================================================================================================#
{#=========================================================#
  #                                                         #
  #                       Combined                          #
  #                                                         #
  #=========================================================#
  data = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/divergence/all_allele_diversity_filt.tab', sep = ';', row.names = NULL)
  data$pi <- data$mpd/data$nseq
  data2 <- merge(data, envir, by = "sample")
  data2$Grouping <- factor(data2$Grouping, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
  data2 <- subset(data2, counts>=10)
  data2$gene <- factor(data2$gene, levels = c("rpoB", "recA", "nodA", "nodD"))
  
  N_data = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/submission_20200918/resubmission_202102/supplementary_data/Table_S1_copy.csv', sep = ";")
  names(N_data)[names(N_data)=="Field.of.origin"] <- "ID"
  
  d <- merge(data2, N_data, by = "ID") # Only DKO matches 
  #============================================#
  #              mpd by effective N            #
  #============================================#
  p <- ggplot( d, aes(x = Effective.N.total..kg.ha., y = mpd, color = gene)) + 
    geom_point() +
    labs( x = "effective N", y = "Pi" ) + 
    theme_bw() 
  p
  
  rpoB <- subset(d, gene=="rpoB")
  correlation(rpoB$mpd, rpoB$Effective.N.total..kg.ha.)
  
  recA <- subset(d, gene=="recA")
  correlation(recA$mpd, recA$Effective.N.total..kg.ha.)
  
  nodA <- subset(d, gene=="nodA")
  correlation(nodA$mpd, nodA$Effective.N.total..kg.ha.)
  
  nodD <- subset(d, gene=="nodD")
  correlation(nodD$mpd, nodD$Effective.N.total..kg.ha.)
}