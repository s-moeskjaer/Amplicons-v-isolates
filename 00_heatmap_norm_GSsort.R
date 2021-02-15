# http://www.roymfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r/

library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)
library(dplyr)
library(viridis)
library(RColorBrewer)


#==========================================#
#                                          #
#        20200304 updated read file        #
#                                          #
#==========================================#

{
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
  
  # Normalize data by row and sample
  
  noda <- nodAfilt[,-22]
  nodd <- nodDfilt[,-18]
  reca <- recAfilt[,-14]
  rpob <- rpoBfilt[,-17]
  
  for (i in 1:nrow(noda)) {
    normalized1 = (noda[i,]/(sum(noda[i,])))
    noda[i,] <- normalized1
  }  
  for (i in 1:nrow(nodd)) {
    normalized2 = (nodd[i,]/(sum(nodd[i,])))
    nodd[i,] <- normalized2
  }  
  for (i in 1:nrow(reca)) {
    normalized3 = (reca[i,]/(sum(reca[i,])))
    reca[i,] <- normalized3
  }  
  for (i in 1:nrow(rpob)) {
    normalized4 = (rpob[i,]/(sum(rpob[i,])))
    rpob[i,] <- normalized4
  }  
  
  # data frame for total gene analysis MERGE
  de <- merge(rpob, reca, by="row.names", all=TRUE)
  rownames(de) <- de$Row.names
  dee <- merge(de, noda, by="row.names", all=TRUE)
  rownames(dee) <- dee$Row.names
  deee <- merge(dee, nodd, by="row.names", all=TRUE)
  rownames(deee) <- deee$Row.names
  deee <- deee[order(as.numeric(deee$Row.names)),]
  
  data_norm <- deee[,-1]
  data_norm <- data_norm[,-1]
  data_norm <- data_norm[,-1]
}



#==========================================#
#                                          #
#               all genes                  #
#                                          #
#==========================================#
{
  #==========================================#
  #                 NA incl                  #
  #==========================================#
  # Reshape from long to wide
  data_norm <- as.data.frame(data_norm)
  data_norm$sampleID <- row.names(data_norm)
  p1.data <- melt(data_norm, id='sampleID')
  
  # Create gene column for facetting 
  p2.data <- data.frame(do.call(rbind, strsplit(as.character(p1.data$variable), "")))
  p3.data <-  paste(p2.data$X1, p2.data$X2, p2.data$X3, p2.data$X4, sep='') 
  p1.data$gene <- p3.data  
  
  # Transform data with log or squareroot 
  p1.data$value <- log10(as.numeric(p1.data$value))
  #p1.data$value <- sqrt(p1.data$value)
  
  
  # Replace NAs with 0
  #p1.data[is.na(p1.data )] <- 0.0001
  #p1.data$value[is.infinite(p1.data$value)] <- 0.0001
  y = subset(p1.data, value > -100 )
  x = min(y$value)
  #p1.data[is.na(p1.data )] <- x
  p1.data$value[is.infinite(p1.data$value)] <- x
  
  # Use pretty row names and discard first column
  row.names(p1.data) <- p1.data$SampleID
  p1.data <- p1.data[ , !(names(p1.data) %in% c('SampleID'))]
  
  # Force order of facets to match list 
  # https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
  p1.data$gene_f = factor(p1.data$gene, levels=c('rpoB','recA','nodA','nodD'))
  p1.data$variable_f = factor(p1.data$variable, levels=c('rpoBseq8', 'rpoBseq9', 'rpoBseq11',
                                                         'rpoBseq13', 'rpoBseq4', 'rpoBseq6', 'rpoBseq10', 'rpoBseq1', 'rpoBseq2', 'rpoBseq3', 'rpoBseq7',
                                                         'rpoBseq12', 'rpoBseq14', 'rpoBseq16', 'rpoBseq5', 'rpoBseq15', 
                                                         'recAseq3', 'recAseq8', 'recAseq2', 'recAseq1', 'recAseq4', 'recAseq6', 'recAseq9', 'recAseq10', 
                                                         'recAseq5', 'recAseq7', 'recAseq11', 'recAseq12', 'recAseq13', 'nodAseq3', 
                                                         'nodAseq10', 'nodAseq2', 'nodAseq6', 'nodAseq11', 'nodAseq12','nodAseq19', 'nodAseq1', 'nodAseq4', 
                                                         'nodAseq8', 'nodAseq16', 'nodAseq22', 'nodAseq9',  'nodAseq5', 'nodAseq13', 'nodAseq15', 'nodAseq18', 
                                                         'nodAseq20', 'nodAseq14', 'nodAseq17', 'nodAseq21', 
                                                         'nodDseq13', 'nodDseq1', 'nodDseq2', 'nodDseq6', 'nodDseq7', 'nodDseq14', 'nodDseq17', 'nodDseq18',
                                                         'nodDseq20', 'nodDseq8', 'nodDseq3','nodDseq5','nodDseq10', 'nodDseq11', 'nodDseq12', 'nodDseq16','nodDseq19'))
  my_palette = brewer.pal(9, "YlGnBu")
  # Plot 
  p = ggplot(p1.data, aes(x=variable_f, y=reorder(sampleID, rev(as.numeric(sampleID))), fill=value)) + 
    geom_tile() +
    #geom_tile(colour="white",size=0.05) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme_bw() +
    facet_grid(~gene_f, switch = "x", scales = "free_x", space = "free_x") + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title=element_text(size = 22),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20), 
          strip.text = element_text(face = "italic")) +
    labs(x="Amplicon", y = "Sample ID") +
    scale_fill_distiller(palette = "YlGnBu", na.value = "gray85") 
  #scale_fill_viridis(na.value = "gray85") # Viridis colour scale
  p
  
  # reverse direction of colour palette: scale_fill_distiller(palette="BrBG", trans = "reverse")
  
  #function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
  ggsave('heatmap_log_facet_NAinc_grey_sort.pdf', plot = p, width = 50, height = 80, unit = 'cm')
  ggsave('heatmap_log_facet_NAinc_grey_sort.png', plot = p, width = 20, height = 30, unit = 'cm')
  
}  



#==========================================#
#                                          #
#            Load all data                 #
#                                          #
#==========================================#

{
  setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/')
  
  data = read.csv('UMI_clover_order.csv', sep = ';')
  head(data)
  names(data)[names(data)=="X..."] <- "sampleID"
  head(data)
  
  #overall sample overview to add field names for DKO
  envir <- read.table('zz_sample_overview.csv', header = T, sep = ";")
  names(envir)[names(envir)=="X...order"] <- "order"
  
  
  ##### FILTER for total UMI count per sample #######
  nodA <- data[,2:23]
  nodD <- data[,24:44]
  recA <- data[,45:57]
  rpoB <- data[,58:73]
  
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
  
  # Normalize data by row and sample
  
  noda <- nodAfilt[,-23]
  nodd <- nodDfilt[,-22]
  reca <- recAfilt[,-14]
  rpob <- rpoBfilt[,-17]
  
  for (i in 1:nrow(noda)) {
    normalized1 = (noda[i,]/(sum(noda[i,])))
    noda[i,] <- normalized1
  }  
  for (i in 1:nrow(nodd)) {
    normalized2 = (nodd[i,]/(sum(nodd[i,])))
    nodd[i,] <- normalized2
  }  
  for (i in 1:nrow(reca)) {
    normalized3 = (reca[i,]/(sum(reca[i,])))
    reca[i,] <- normalized3
  }  
  for (i in 1:nrow(rpob)) {
    normalized4 = (rpob[i,]/(sum(rpob[i,])))
    rpob[i,] <- normalized4
  }  
  
  # data frame for total gene analysis MERGE
  de <- merge(rpob, reca, by="row.names", all=TRUE)
  rownames(de) <- de$Row.names
  dee <- merge(de, noda, by="row.names", all=TRUE)
  rownames(dee) <- dee$Row.names
  deee <- merge(dee, nodd, by="row.names", all=TRUE)
  rownames(deee) <- deee$Row.names
  deee <- deee[order(as.numeric(deee$Row.names)),]
  
  data_norm <- deee[,-1]
  data_norm <- data_norm[,-1]
  data_norm <- data_norm[,-1]
}



#==========================================#
#                                          #
#               all genes                  #
#                                          #
#==========================================#
{
  #==========================================#
  #                 NA incl                  #
  #==========================================#
  # Reshape from long to wide
  data_norm <- as.data.frame(data_norm)
  data_norm$sampleID <- row.names(data_norm)
  p1.data <- melt(data_norm, id='sampleID')
  
  # Create gene column for facetting 
  p2.data <- data.frame(do.call(rbind, strsplit(as.character(p1.data$variable), "")))
  p3.data <-  paste(p2.data$X1, p2.data$X2, p2.data$X3, p2.data$X4, sep='') 
  p1.data$gene <- p3.data  
  
  # Transform data with log or squareroot 
  p1.data$value <- log10(as.numeric(p1.data$value))
  #p1.data$value <- sqrt(p1.data$value)
  
  
  # Replace NAs with 0
  #p1.data[is.na(p1.data )] <- 0.0001
  #p1.data$value[is.infinite(p1.data$value)] <- 0.0001
  y = subset(p1.data, value > -100 )
  x = min(y$value)
  p1.data[is.na(p1.data )] <- x
  p1.data$value[is.infinite(p1.data$value)] <- x
  
  # Use pretty row names and discard first column
  row.names(p1.data) <- p1.data$SampleID
  p1.data <- p1.data[ , !(names(p1.data) %in% c('SampleID'))]
  
  # Force order of facets to match list 
  # https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
  p1.data$gene_f = factor(p1.data$gene, levels=c('rpoB','recA','nodA','nodD'))
  
  my_palette = brewer.pal(9, "YlGnBu")
  # Plot 
  p = ggplot(p1.data, aes(x=variable, y=reorder(sampleID, rev(as.numeric(sampleID))), fill=value)) + 
    geom_tile() +
    #geom_tile(colour="white",size=0.05) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme_bw() +
    facet_grid(~gene_f, switch = "x", scales = "free_x", space = "free_x") + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title=element_text(size = 22),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20), 
          strip.text = element_text(face = "italic")) +
    labs(x="Amplicon", y = "Sample ID") +
    scale_fill_distiller(palette = "YlGnBu") 
  #scale_fill_viridis(na.value = "gray85") # Viridis colour scale
  p
  
  # reverse direction of colour palette: scale_fill_distiller(palette="BrBG", trans = "reverse")
  
  #function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
  ggsave('heatmap_log_facet_NAinc_blue.pdf', plot = p, width = 50, height = 80, unit = 'cm')
  
  
  
  
  #==========================================#
  #                 NA omit                  #
  #==========================================#
  data_norm <- deee[,-1]
  data_norm <- data_norm[,-1]
  data_norm <- data_norm[,-1]
  data_norm <- na.omit(data_norm)
  # Reshape from long to wide
  data_norm <- as.data.frame(data_norm)
  data_norm$sampleID <- row.names(data_norm)
  p1.data <- melt(data_norm, id='sampleID')
  
  # Create gene column for facetting 
  p2.data <- data.frame(do.call(rbind, strsplit(as.character(p1.data$variable), "")))
  p3.data <-  paste(p2.data$X1, p2.data$X2, p2.data$X3, p2.data$X4, sep='') 
  p1.data$gene <- p3.data  
  
  # Transform data with log or squareroot 
  p1.data$value <- log10(as.numeric(p1.data$value))
  #p1.data$value <- sqrt(p1.data$value)
  
  
  # Replace NAs with 0
  #p1.data[is.na(p1.data )] <- 0.0001
  #p1.data$value[is.infinite(p1.data$value)] <- 0.0001
  y = subset(p1.data, value > -100 )
  x = min(y$value)
  #p1.data[is.na(p1.data )] <- x
  p1.data$value[is.infinite(p1.data$value)] <- x
  
  # Use pretty row names and discard first column
  row.names(p1.data) <- p1.data$SampleID
  p1.data <- p1.data[ , !(names(p1.data) %in% c('SampleID'))]
  
  # Force order of facets to match list 
  # https://stackoverflow.com/questions/14262497/fixing-the-order-of-facets-in-ggplot
  p1.data$gene_f = factor(p1.data$gene, levels=c('rpoB','recA','nodA','nodD'))
  
  # Plot 
  p = ggplot(p1.data, aes(x=variable, y=reorder(sampleID, rev(as.numeric(sampleID))), fill=value)) + 
    geom_tile() +
    #geom_tile(colour="white",size=0.05) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme_bw() +
    facet_grid(~gene_f, switch = "x", scales = "free_x", space = "free_x") + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title=element_text(size = 22),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20), 
          strip.text = element_text(face = "italic")) +
    labs(x="Amplicon", y = "Sample ID") +
    scale_fill_viridis()
  p
  #function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
  ggsave('heatmap_log_facet_NAomit.pdf', plot = p, width = 50, height = 80, unit = 'cm')
}




#==========================================#
#                                          #
#                   recA                   #
#                                          #
#==========================================#

{
  data_norm <- reca
  
  # Add sampleID back in data frame 
  data_norm$sampleID <- row.names(data_norm)
  data_norm$index <- as.numeric(row.names(data_norm))
  data_norm <- data_norm[order(data_norm$index), ]
  data_norm$index <- NULL
  # Reshape from long to wide
  data_norm <- as.data.frame(data_norm)
  #data <- data[1:32,] 
  p1.data <- melt(data_norm, id='sampleID')
  
  # Transform data with log or squareroot 
  p1.data$value <- log10(as.numeric(p1.data$value))
  #p1.data$value <- sqrt(p1.data$value)
  
  
  # Replace -infinite values with lowest value in dataframe
  y = subset(p1.data, value > -100 )
  x = min(y$value)
  p1.data$value[is.infinite(p1.data$value)] <- x
  
  # Use pretty row names and discard first column
  row.names(p1.data) <- p1.data$SampleID
  p1.data <- p1.data[ , !(names(p1.data) %in% c('SampleID'))]
  
  # Plot 
  p = ggplot(p1.data, aes(x=variable, y=reorder(sampleID, rev(as.numeric(sampleID))), fill=value)) + 
    geom_tile() +
    #geom_tile(colour="white",size=0.05) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title=element_text(size = 22),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20), 
          strip.text = element_text(face = "italic")) +
    labs(x="Amplicon", y = "Sample ID") +
    scale_fill_viridis()
  p
  
  #additional parameters(width=,height=,units='cm')
  ggsave('heatmap_log_recA.pdf', plot = p, width = ncol(data_norm), height = 50, unit = 'cm')
}

#==========================================#
#                                          #
#                   rpob                   #
#                                          #
#==========================================#

{
  data_norm <- rpob
  
  # Add sampleID back in data frame 
  data_norm$sampleID <- row.names(data_norm)
  data_norm$index <- as.numeric(row.names(data_norm))
  data_norm <- data_norm[order(data_norm$index), ]
  data_norm$index <- NULL
  # Reshape from long to wide
  data_norm <- as.data.frame(data_norm)
  #data <- data[1:32,] 
  p1.data <- melt(data_norm, id='sampleID')
  
  # Transform data with log or squareroot 
  p1.data$value <- log10(as.numeric(p1.data$value))
  #p1.data$value <- sqrt(p1.data$value)
  
  
  # Replace -infinite values with lowest value in dataframe
  y = subset(p1.data, value > -100 )
  x = min(y$value)
  p1.data$value[is.infinite(p1.data$value)] <- x
  
  # Use pretty row names and discard first column
  row.names(p1.data) <- p1.data$SampleID
  p1.data <- p1.data[ , !(names(p1.data) %in% c('SampleID'))]
  
  # Plot 
  p = ggplot(p1.data, aes(x=variable, y=reorder(sampleID, rev(as.numeric(sampleID))), fill=value)) + 
    geom_tile() +
    #geom_tile(colour="white",size=0.05) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title=element_text(size = 22),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20), 
          strip.text = element_text(face = "italic")) +
    labs(x="Amplicon", y = "Sample ID") +
    scale_fill_viridis()
  p
  #function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
  ggsave('heatmap_log_rpob.pdf', plot = p, width = ncol(data_norm), height = 50, unit = 'cm')
}



#==========================================#
#                                          #
#                   nodA                   #
#                                          #
#==========================================#


{
  data_norm <- noda
  
  # Add sampleID back in data frame 
  data_norm$sampleID <- row.names(data_norm)
  data_norm$index <- as.numeric(row.names(data_norm))
  data_norm <- data_norm[order(data_norm$index), ]
  data_norm$index <- NULL
  # Reshape from long to wide
  data_norm <- as.data.frame(data_norm)
  #data <- data[1:32,] 
  p1.data <- melt(data_norm, id='sampleID')
  
  # Transform data with log or squareroot 
  p1.data$value <- log10(as.numeric(p1.data$value))
  #p1.data$value <- sqrt(p1.data$value)
  
  
  # Replace -infinite values with lowest value in dataframe
  y = subset(p1.data, value > -100 )
  x = min(y$value)
  p1.data$value[is.infinite(p1.data$value)] <- x
  
  # Use pretty row names and discard first column
  row.names(p1.data) <- p1.data$SampleID
  p1.data <- p1.data[ , !(names(p1.data) %in% c('SampleID'))]
  
  # Plot 
  p = ggplot(p1.data, aes(x=variable, y=reorder(sampleID, rev(as.numeric(sampleID))), fill=value)) + 
    geom_tile() +
    #geom_tile(colour="white",size=0.05) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title=element_text(size = 22),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20), 
          strip.text = element_text(face = "italic")) +
    labs(x="Amplicon", y = "Sample ID") +
    scale_fill_viridis()
  p
  #function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
  ggsave('heatmap_log_nodA.pdf', plot = p, width = ncol(data_norm), height = 80, unit = 'cm')
}



#==========================================#
#                                          #
#                   nodD                   #
#                                          #
#==========================================#
{
  data_norm <- nodd
  
  # Add sampleID back in data frame 
  data_norm$sampleID <- row.names(data_norm)
  data_norm$index <- as.numeric(row.names(data_norm))
  data_norm <- data_norm[order(data_norm$index), ]
  data_norm$index <- NULL
  # Reshape from long to wide
  data_norm <- as.data.frame(data_norm)
  #data <- data[1:32,] 
  p1.data <- melt(data_norm, id='sampleID')
  
  # Transform data with log or squareroot 
  p1.data$value <- log10(as.numeric(p1.data$value))
  #p1.data$value <- sqrt(p1.data$value)
  
  
  # Replace -infinite values with lowest value in dataframe
  y = subset(p1.data, value > -100 )
  x = min(y$value)
  p1.data$value[is.infinite(p1.data$value)] <- x
  
  # Use pretty row names and discard first column
  row.names(p1.data) <- p1.data$SampleID
  p1.data <- p1.data[ , !(names(p1.data) %in% c('SampleID'))]
  
  # Plot 
  p = ggplot(p1.data, aes(x=variable, y=reorder(sampleID, rev(as.numeric(sampleID))), fill=value)) + 
    geom_tile() +
    #geom_tile(colour="white",size=0.05) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title=element_text(size = 22),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 20), 
          strip.text = element_text(face = "italic")) +
    labs(x="Amplicon", y = "Sample ID") +
    scale_fill_viridis()
  p
  #function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
  ggsave('heatmap_log_nodd.pdf', plot = p, width = ncol(data_norm), height = 80, unit = 'cm')
}
