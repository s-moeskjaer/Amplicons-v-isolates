setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/divergence/')
library(ggplot2)
library(dplyr)
library(multcompView)


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
  
  #=============================#
  #              mpd            #
  #=============================#
  plot2 <- ggplot(data2, aes(x=gene, y=pi, fill = Grouping)) +
    stat_boxplot(geom = "errorbar", position=position_dodge(0.9)) +
    geom_boxplot(outlier.size=0.6,color = 'black', position=position_dodge(0.9)) +
    geom_point(position=position_dodge(0.9), size=0.6) +
    theme_bw() +
    labs(x="Amplicon", y = ~ pi) +
    #scale_fill_manual(values=c("#b2182b", "#f4a582", "#e0e0e0", "#4d4d4d")) # grey-red colours
    #scale_fill_manual(values=c("#0f1027", "#5c2e88", "#ba4f8f", "#f68f75")) + # petroleum-brown colours
    scale_fill_manual(values=c("#40004b", "#762a83", "#c2a5cf", "#b8e186", "#4d9221", "#276419", "#fed976", "#fd8d3c", "#bd0026")) +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(face="italic", size = 20),
          axis.text.y = element_text(size = 20),
          axis.title=element_text(size = 22),
          legend.text = element_text(size = 20))
  
  plot2
  ggsave('202007_all_mpd_filt_dots.pdf', plot = plot2, width = 30, height = 20, unit = 'cm', useDingbats = FALSE)
  #ggsave('202007_all_mpd_filt.png', plot = plot2, width = 30, height = 20, unit = 'cm')
  
  #=========================================================#
  #                                                         #
  #                      DLF v DKO                          #
  #                                                         #
  #=========================================================#
  envir <- read.table('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/divergence/divergence_DKOvDLF/zz_sample_overview.csv', header = T, sep = ";")

  data2 <- merge(data, envir, by = "sample")
  data2$Management <- factor(data2$Management, levels = c("DKO", "DLF"))
  data2 <- subset(data2, counts>=10)
  data2$gene <- factor(data2$gene, levels = c("rpoB", "recA", "nodA", "nodD"))
  
  #=============================#
  #              mpd            #
  #=============================#
  plot2 <- ggplot(data2, aes(x=gene, y=pi, fill = Management)) +
    stat_boxplot(geom = "errorbar") +
    geom_boxplot(outlier.size=0.5,color = 'black') +
    geom_point(size=0.6, position=position_dodge(0.75)) +
    theme_bw() +
    labs(x="Amplicon", y = ~ pi) +
    #scale_fill_manual(values=c("#b2182b", "#f4a582", "#e0e0e0", "#4d4d4d")) # grey-red colours
    #scale_fill_manual(values=c("#0f1027", "#5c2e88", "#ba4f8f", "#f68f75")) + # petroleum-brown colours
    scale_fill_manual(values=c("#d6604d", "#4d4d4d")) +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(face="italic", size = 20),
          axis.text.y = element_text(size = 20),
          axis.title=element_text(size = 22),
          legend.text = element_text(size = 20))
  
  
  plot2
  ggsave('202007_all_mpd_filt_DLFvDKO_dots.pdf', plot = plot2, width = 30, height = 20, unit = 'cm', useDingbats = FALSE)
  
  # ANOVA and t test for DLF v DKO
  nodD <- subset(data2, gene == "nodD")
  fit <- aov(pi ~ Management  , data=nodD)
  results <- TukeyHSD(fit)
  multcompLetters4(fit, results)
  t.test(nodD$pi~nodD$Management)
  
  
  nodA <- subset(data2, gene == "nodA")
  fit <- aov(pi ~ Management  , data=nodA)
  results <- TukeyHSD(fit)
  multcompLetters4(fit, results)
  t.test(nodA$pi~nodA$Management)
  
  rpoB <- subset(data2, gene == "rpoB")
  fit <- aov(pi ~ Management  , data=rpoB)
  results <- TukeyHSD(fit)
  multcompLetters4(fit, results)
  t.test(rpoB$pi~rpoB$Management)
  
  recA <- subset(data2, gene == "recA")
  fit <- aov(pi ~ Management  , data=recA)
  results <- TukeyHSD(fit)
  multcompLetters4(fit, results)
  t.test(recA$pi~recA$Management)
  
  
  nrow(subset(nodD, Management == "DKO"))
  nrow(subset(nodD, Management == "DLF"))
  
  # ANOVA and t test for Groupings
  nodD <- subset(data2, gene == "nodD")
  fit <- aov(pi ~ Grouping  , data=nodD)
  results <- TukeyHSD(fit)
  multcompLetters4(fit, results)
  
  nodA <- subset(data2, gene == "nodA")
  fit <- aov(pi ~ Grouping  , data=nodA)
  results <- TukeyHSD(fit)
  multcompLetters4(fit, results)
  
  rpoB <- subset(data2, gene == "rpoB")
  fit <- aov(pi ~ Grouping  , data=rpoB)
  results <- TukeyHSD(fit)
  multcompLetters4(fit, results)
  
  
  recA <- subset(data2, gene == "recA")
  fit <- aov(pi ~ Grouping  , data=recA)
  results <- TukeyHSD(fit)
  multcompLetters4(fit, results)
}


