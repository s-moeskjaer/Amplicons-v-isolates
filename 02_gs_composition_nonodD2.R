#==========================================#
#                                          #
#        Initial data normalization        #
#                                          #
#==========================================#
{
  library(reshape2)
  library(ggplot2)
  library("lattice")
  setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/')
  data = read.csv('UMI_clover_order_nonodD2.csv', sep = ';')
  head(data)
  names(data)[names(data)=="X"] <- "sampleID"
  head(data)
  
  #overall sample overview to add field names for DKO
  envir <- read.table('zz_sample_overview.csv', header = T, sep = ";")

  
  ##### FILTER for total UMI count per sample #######
  nodA <- data[,2:22]
  nodD <- data[,23:39]
  recA <- data[,40:52]
  rpoB <- data[,53:68]
  
  nodA$total <- rowSums(nodA)
  rownames(nodA) <- data$sampleID
  nodD$total <- rowSums(nodD)
  rownames(nodD) <- data$sampleID
  recA$total <- rowSums(recA)
  rownames(recA) <- data$sampleID
  rpoB$total <- rowSums(rpoB)
  rownames(rpoB) <- data$sampleID
  
  # Filter out samples with total UMI count<10 
  nodAfilt <- subset(nodA, total>=10)
  nodDfilt <- subset(nodD, total>=10)
  recAfilt <- subset(recA, total>=10)
  rpoBfilt <- subset(rpoB, total>=10)
  
  nodAfilt <- nodAfilt[,-22]
  nodDfilt <- nodDfilt[,-18]
  recAfilt <- recAfilt[,-14]
  rpoBfilt <- rpoBfilt[,-17]
  
}


#==========================================#
#                                          #
#                   rpoB                   #
#                                          #
#==========================================#
#==========================================#
#             DKO groupings                #
#==========================================#
{
genoPalette <- c("#466EA9","#EA877A", "#4D9A7A", "#8E569D", "#DA367E")

rpob_samples <- envir
row.names(rpob_samples) <- rpob_samples$order
d1 <- merge(rpoBfilt,rpob_samples, by="row.names",all.x=FALSE)
d2 <- d1[,-1]
d2 <- d2[,-19:-23]
d2 <- d2[,-17]

DF1 <- melt(d2, id.var=c("Sample", "Grouping"))

#Add genospecies column 
genospecies = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/amplicon_strainID_nonodD2.csv', sep = ';')
names(genospecies)[names(genospecies)=="amplicon"] <- "variable"
gs <- genospecies[,-2:-11]
DF2 <- merge(DF1, gs, by = "variable")
# Force order of facets to match list 
# https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
DF2$Grouping2 = factor(DF2$Grouping, levels=c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))

# scaled to percentage 
p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
  expand_limits(x = 0, y = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  facet_grid(~Grouping2, switch = "x", scales = "free_x", space = "free_x") + 
  labs(y = "Genospecies abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 22),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20))

ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_rpoB.pdf', plot = p, width = 30, height = 10, unit = 'cm')


# true numbers 
p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
  expand_limits(x = 0, y = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  facet_grid(~Grouping2, switch = "x", scales = "free_x", space = "free_x") + 
  labs(y = "Genospecies abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 22),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20))

p

ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_raw_rpoB.pdf', plot = p, width = 30, height = 10, unit = 'cm')
}
#==========================================#
#                 DKO DLF                  #
#==========================================#
{
genoPalette <- c("#466EA9","#EA877A", "#4D9A7A", "#8E569D", "#DA367E")

rpob_samples <- envir
row.names(rpob_samples) <- rpob_samples$order
d1 <- merge(rpoBfilt,rpob_samples, by="row.names",all.x=FALSE)
d2 <- d1[,-1]
d2 <- d2[,-21:-23]
d2 <- d2[,-19]
d2 <- d2[,-17]

DF1 <- melt(d2, id.var=c("Sample", "Origin", "Grouping"))


#Add genospecies column 
genospecies = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/amplicon_strainID_nonodD2.csv', sep = ';')
names(genospecies)[names(genospecies)=="amplicon"] <- "variable"
gs <- genospecies[,-2:-11]
DF2 <- merge(DF1, gs, by = "variable")
# Force order of facets to match list 
# https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
DF2$Origin = factor(DF2$Origin, levels=c("DKO", "DK", "F", "UK"))

# scaled to percentage 
p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
  expand_limits(x = 0, y = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  facet_grid(~Origin, switch = "x", scales = "free_x", space = "free_x") + 
  labs(y = "Genospecies abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 22),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20))

ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_rpoB_DKO_DLF.pdf', plot = p, width = 30, height = 10, unit = 'cm')


# true numbers 
p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
  expand_limits(x = 0, y = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  facet_grid(~Origin, switch = "x", scales = "free_x", space = "free_x") + 
  labs(y = "Genospecies abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 22),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20))

ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_raw_rpoB_DKO_DLF.pdf', plot = p, width = 30, height = 10, unit = 'cm')
}


#==========================================#
#                                          #
#                   recA                   #
#                                          #
#==========================================#

#==========================================#
#             DKO groupings                #
#==========================================#
{
genoPalette <- c("#466EA9","#EA877A", "#4D9A7A", "#8E569D", "#DA367E")

reca_samples <- envir
row.names(reca_samples) <- reca_samples$order
d1 <- merge(recAfilt,reca_samples, by="row.names",all.x=FALSE)
d2 <- d1[,-1]
d2 <- d2[,-16:-20]
d2 <- d2[,-14]

DF1 <- melt(d2, id.var=c("Sample", "Grouping"))

#Add genospecies column 
genospecies = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/amplicon_strainID_nonodD2.csv', sep = ';')
names(genospecies)[names(genospecies)=="amplicon"] <- "variable"
gs <- genospecies[,-2:-11]
DF2 <- merge(DF1, gs, by = "variable")
# Force order of facets to match list 
# https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
DF2$Grouping2 = factor(DF2$Grouping, levels=c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))

# scaled to percentage 
p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
  expand_limits(x = 0, y = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  facet_grid(~Grouping2, switch = "x", scales = "free_x", space = "free_x") + 
  labs(y = "Genospecies abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 22),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20))

ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_recA.pdf', plot = p, width = 30, height = 10, unit = 'cm')


# true numbers 
p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
  expand_limits(x = 0, y = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  facet_grid(~Grouping2, switch = "x", scales = "free_x", space = "free_x") + 
  labs(y = "Genospecies abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 22),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20))

ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_raw_recA.pdf', plot = p, width = 30, height = 10, unit = 'cm')
}

#==========================================#
#                 DKO DLF                  #
#==========================================#
{
genoPalette <- c("#466EA9","#EA877A", "#4D9A7A", "#8E569D", "#DA367E")

reca_samples <- envir
row.names(reca_samples) <- reca_samples$order
d1 <- merge(recAfilt,reca_samples, by="row.names",all.x=FALSE)
d2 <- d1[,-1]
d2 <- d2[,-18:-21]
d2 <- d2[,-16]
d2 <- d2[,-14]

DF1 <- melt(d2, id.var=c("Sample", "Origin"))

#Add genospecies column 
genospecies = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/amplicon_strainID_nonodD2.csv', sep = ';')
names(genospecies)[names(genospecies)=="X...amplicon"] <- "variable"
gs <- genospecies[,-2:-11]
DF2 <- merge(DF1, gs, by = "variable")
# Force order of facets to match list 
# https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
DF2$Origin = factor(DF2$Origin, levels=c("DKO", "DK", "F", "UK"))

# scaled to percentage 
p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
  expand_limits(x = 0, y = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  facet_grid(~Origin, switch = "x", scales = "free_x", space = "free_x") + 
  labs(y = "Genospecies abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 22),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20))

ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_recA_DKO_DLF.pdf', plot = p, width = 30, height = 10, unit = 'cm')


# true numbers 
p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
  expand_limits(x = 0, y = 0) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  facet_grid(~Origin, switch = "x", scales = "free_x", space = "free_x") + 
  labs(y = "Genospecies abundance") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 22),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20))

ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_raw_recA_DKO_DLF.pdf', plot = p, width = 30, height = 10, unit = 'cm')

}


#==========================================#
#                                          #
#                   nodA                   #
#                                          #
#==========================================#
#==========================================#
#             DKO groupings                #
#==========================================#
{
  genoPalette <- c("#466EA9","#EA877A", "#4D9A7A", "#DA367E", "#e0f3db")
  
  noda_samples <- envir
  row.names(noda_samples) <- noda_samples$order
  d1 <- merge(nodAfilt,noda_samples, by="row.names",all.x=FALSE)
  d2 <- d1[,-1]
  d2 <- d2[,-24:-28]
  d2 <- d2[,-22]
  
  DF1 <- melt(d2, id.var=c("Sample", "Grouping"))

  #Add genospecies column 
  genospecies = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/amplicon_strainID_nonodD2.csv', sep = ';')
  names(genospecies)[names(genospecies)=="amplicon"] <- "variable"
  gs <- genospecies[,-2:-11]
  DF2 <- merge(DF1, gs, by = "variable")
  # Force order of facets to match list 
  # https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
  DF2$Grouping2 = factor(DF2$Grouping, levels=c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
  
  # scaled to percentage 
  p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    facet_grid(~Grouping2, switch = "x", scales = "free_x", space = "free_x") + 
    labs(y = "Genospecies abundance") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size = 22),
          axis.title.x=element_blank(),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20))
  
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_nodA.pdf', plot = p, width = 30, height = 10, unit = 'cm')
  
  
  # true numbers 
  p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) + 
    facet_grid(~Grouping2, switch = "x", scales = "free_x", space = "free_x") + 
    labs(y = "Genospecies abundance") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size = 22),
          axis.title.x=element_blank(),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20))
  
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_raw_nodA.pdf', plot = p, width = 30, height = 10, unit = 'cm')
}

#==========================================#
#                 DKO DLF                  #
#==========================================#
{
  genoPalette <- c("#466EA9","#EA877A", "#4D9A7A", "#DA367E", "#e0f3db")
  
  noda_samples <- envir
  row.names(noda_samples) <- noda_samples$order
  d1 <- merge(nodAfilt,noda_samples, by="row.names",all.x=FALSE)
  d2 <- d1[,-1]
  d2 <- d2[,-26:-29]
  d2 <- d2[,-24]
  d2 <- d2[,-22]
  
  DF1 <- melt(d2, id.var=c("Sample", "Origin"))
  
  #Add genospecies column 
  genospecies = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/amplicon_strainID_nonodD2.csv', sep = ';')
  names(genospecies)[names(genospecies)=="X...amplicon"] <- "variable"
  gs <- genospecies[,-2:-11]
  DF2 <- merge(DF1, gs, by = "variable")
  # Force order of facets to match list 
  # https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
  DF2$Origin = factor(DF2$Origin, levels=c("DKO", "DK", "F", "UK"))
  
  # scaled to percentage 
  p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    facet_grid(~Origin, switch = "x", scales = "free_x", space = "free_x") + 
    labs(y = "Genospecies abundance") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size = 22),
          axis.title.x=element_blank(),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20))
  
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_nodA_DKO_DLF.pdf', plot = p, width = 30, height = 10, unit = 'cm')
  
  
  # true numbers 
  p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    facet_grid(~Origin, switch = "x", scales = "free_x", space = "free_x") + 
    labs(y = "Genospecies abundance") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size = 22),
          axis.title.x=element_blank(),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20))
  
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_raw_nodA_DKO_DLF.pdf', plot = p, width = 30, height = 10, unit = 'cm')
  
}

#==========================================#
#                                          #
#                   nodD                   #
#                                          #
#==========================================#
#==========================================#
#             DKO groupings                #
#==========================================#
{
  genoPalette <- c("#466EA9","#EA877A", "#4D9A7A", "#DA367E", "#e0f3db")
  
  nodd_samples <- envir
  row.names(nodd_samples) <- nodd_samples$order
  d1 <- merge(nodDfilt,nodd_samples, by="row.names",all.x=FALSE)
  d2 <- d1[,-1]
  d2 <- d2[,-20:-24]
  d2 <- d2[,-18]
  
  DF1 <- melt(d2, id.var=c("Sample", "Grouping"))
  
  #Add genospecies column 
  genospecies = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/amplicon_strainID_nonodD2.csv', sep = ';')
  names(genospecies)[names(genospecies)=="amplicon"] <- "variable"
  gs <- genospecies[,-2:-11]
  DF2 <- merge(DF1, gs, by = "variable")
  # Force order of facets to match list 
  # https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
  DF2$Grouping2 = factor(DF2$Grouping, levels=c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
  
  # scaled to percentage 
  p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    facet_grid(~Grouping2, switch = "x", scales = "free_x", space = "free_x") + 
    labs(y = "Genospecies abundance") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size = 22),
          axis.title.x=element_blank(),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20))
  
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_nodD.pdf', plot = p, width = 30, height = 10, unit = 'cm')
  
  
  # true numbers 
  p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) + 
    facet_grid(~Grouping2, switch = "x", scales = "free_x", space = "free_x") + 
    labs(y = "Genospecies abundance") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size = 22),
          axis.title.x=element_blank(),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20))
  
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_raw_nodD.pdf', plot = p, width = 30, height = 10, unit = 'cm')
}

#==========================================#
#                 DKO DLF                  #
#==========================================#
{
  genoPalette <- c("#466EA9","#EA877A", "#4D9A7A", "#DA367E", "#e0f3db")
  
  nodd_samples <- envir
  row.names(nodd_samples) <- nodd_samples$order
  d1 <- merge(nodDfilt,nodd_samples, by="row.names",all.x=FALSE)
  d2 <- d1[,-1]
  d2 <- d2[,-22:-25]
  d2 <- d2[,-18]
  d2 <- d2[,-19]
  
  DF1 <- melt(d2, id.var=c("Sample", "Origin"))
  
  #Add genospecies column 
  genospecies = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/amplicon_strainID_nonodD2.csv', sep = ';')
  names(genospecies)[names(genospecies)=="X...amplicon"] <- "variable"
  gs <- genospecies[,-2:-11]
  DF2 <- merge(DF1, gs, by = "variable")
  # Force order of facets to match list 
  # https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
  DF2$Origin = factor(DF2$Origin, levels=c("DKO", "DK", "F", "UK"))
  
  # scaled to percentage 
  p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    facet_grid(~Origin, switch = "x", scales = "free_x", space = "free_x") + 
    labs(y = "Genospecies abundance") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size = 22),
          axis.title.x=element_blank(),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20))
  
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_nodD_DKO_DLF.pdf', plot = p, width = 30, height = 10, unit = 'cm')
  
  
  # true numbers 
  p = ggplot(DF2, aes(x = Sample, y = value, fill = phylogeny)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual("legend", values = genoPalette,  na.value = "#6d6e71") +
    expand_limits(x = 0, y = 0) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    facet_grid(~Origin, switch = "x", scales = "free_x", space = "free_x") + 
    labs(y = "Genospecies abundance") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size = 22),
          axis.title.x=element_blank(),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20))
  
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/gs_comp/no_nodD2_phyl/gs_comp_raw_nodD_DKO_DLF.pdf', plot = p, width = 30, height = 10, unit = 'cm')
  
}

