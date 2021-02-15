
library(ggplot2)  
library(reshape2)
library(plyr)
library(data.table)
library(dplyr)
library(viridis)
library(BBmisc)
library(caret)
library(ggfortify)
library(ggpubr)
 

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



#=======================================================#
#                                                       #
#                   By Grouping                         #
#                                                       #
#=======================================================#
#==========================================#
#                ALL NA incl               #
#==========================================#
{
  data.pca <- data_norm
  # Replace 0s and NA with 0.0000000001
  data.pca[data.pca == 0] <- 0.0000000001
  data.pca[is.na(data.pca )] <- 0.0000000001

  samples <- row.names(data_norm)
  
  #Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  
  # Create plot column for colouring
  plot.names <- d2$Grouping
  
  # Log transform data for pca
  # Change end row value to be the ncolumns
  x = ncol(data.pca)
  data.pca[ ,1:x] <- log10(data.pca[ ,1:x])
  
  # Remove 0 variance using caret for when Aarhus plots are removed previously
  NZV <- nearZeroVar(data.pca, saveMetrics = TRUE)
  nzv.col <- nearZeroVar(data.pca)
  if(length(nzv.col) > 0) data.pca <- data.pca[, -nzv.col]
  
  d.pca <- prcomp(data.pca,
                  center = TRUE,
                  scale. = TRUE)
  
  # Create percent variance label for axis 
  percentVar <- d.pca$sdev^2 / sum( d.pca$sdev^2 )
  percentVar <- round(100 * percentVar)
  
  #Order groupings 
  plot.names <- factor(plot.names, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
  grouping.colours <- c("#40004b", "#762a83", "#c2a5cf", "#b8e186", "#4d9221", "#276419", "#fed976", "#fd8d3c", "#bd0026")
  
  (pcaplot <- ggplot(d.pca, aes(PC1, PC2, color = plot.names , fill = plot.names)) +
      coord_fixed() +
      geom_point(size=3) +
      scale_fill_manual(values = grouping.colours) +
      scale_colour_manual(values = grouping.colours) +
      #must manually apply the percentage variance in this case
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      #geom_polygon(data = hulls, alpha = 0.5) +
      stat_chull(geom = "polygon", alpha = 0.2) +
      #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
      #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
      stat_mean(size = 5) +
      theme_bw() +
      theme(aspect.ratio=1,
            legend.title = element_blank(),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title=element_text(size = 22),
            legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_all_PC12.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm' , useDingbats=FALSE)
  
  ########################
  ######### PC23 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC2, PC3, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC2: ",percentVar[2],"% variance")) +
     ylab(paste0("PC3: ",percentVar[3],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_all_PC23.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
  
  ########################
  ######### PC34 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC3, PC4, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #Manually apply the percentage variance in this case
     xlab(paste0("PC3: ",percentVar[3],"% variance")) +
     ylab(paste0("PC4: ",percentVar[4],"% variance")) + 
     stat_chull(geom = "polygon", alpha = 0.2) +
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_all_PC34.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  

}
#==========================================#
#                ALL NA omit               #
#==========================================#
{
  data.pca <- na.omit(data_norm)
  # Replace 0s and NA with 0.0000000001
  data.pca[data.pca == 0] <- 0.0000000001

  samples <- row.names(data.pca)
  
  #Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  
  # Create plot column for colouring
  plot.names <- d2$Grouping
  
  # Log transform data for pca
  # Change end row value to be the ncolumns
  x = ncol(data.pca)
  data.pca[ ,1:x] <- log10(data.pca[ ,1:x])
  
  # Remove 0 variance using caret for when Aarhus plots are removed previously
  NZV <- nearZeroVar(data.pca, saveMetrics = TRUE)
  nzv.col <- nearZeroVar(data.pca)
  if(length(nzv.col) > 0) data.pca <- data.pca[, -nzv.col]
  
  d.pca <- prcomp(data.pca,
                  center = TRUE,
                  scale. = TRUE)
  
  # Create percent variance label for axis 
  percentVar <- d.pca$sdev^2 / sum( d.pca$sdev^2 )
  percentVar <- round(100 * percentVar)
  
  #Order groupings 
  plot.names <- factor(plot.names, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
  grouping.colours <- c("#40004b", "#762a83", "#c2a5cf", "#b8e186", "#4d9221", "#276419", "#fed976", "#fd8d3c", "#bd0026")
  
  (pcaplot <- ggplot(d.pca, aes(PC1, PC2, color = plot.names , fill = plot.names)) +
      coord_fixed() +
      geom_point(size=3) +
      scale_fill_manual(values = grouping.colours) +
      scale_colour_manual(values = grouping.colours) +
      #must manually apply the percentage variance in this case
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      #geom_polygon(data = hulls, alpha = 0.5) +
      stat_chull(geom = "polygon", alpha = 0.2) +
      #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
      #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
      stat_mean(size = 5) +
      theme_bw() +
      theme(aspect.ratio=1,
            legend.title = element_blank(),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title=element_text(size = 22),
            legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_all_naomit_PC12.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm' , useDingbats=FALSE)
  
  ########################
  ######### PC23 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC2, PC3, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC2: ",percentVar[2],"% variance")) +
     ylab(paste0("PC3: ",percentVar[3],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_all_naomit_PC23.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
  
  ########################
  ######### PC34 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC3, PC4, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC3: ",percentVar[3],"% variance")) +
     ylab(paste0("PC4: ",percentVar[4],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_all_naomit_PC34.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
}



#================================================#
#                    rpoB                        #
#================================================#
{
  #Order 
  rpob$index <- as.numeric(row.names(rpob))
  rpob <- rpob[order(rpob$index), ]
  samples <- rpob$index
  rpob$index <- NULL
  
  #Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  
  # Create plot column for colouring
  plot.names <- d2$Grouping
  
  data.pca <- rpob
  
  # Replace 0s with 0.00001
  data.pca[data.pca == 0] <- 0.0000000001
  
  # Log transform data for pca
  # Change end row value to be the ncolumns
  x = ncol(data.pca)
  data.pca[ ,1:x] <- log10(data.pca[ ,1:x])
  
  # Remove 0 variance using caret for when Aarhus plots are removed previously
  NZV <- nearZeroVar(data.pca, saveMetrics = TRUE)
  nzv.col <- nearZeroVar(data.pca)
  if(length(nzv.col) > 0) data.pca <- data.pca[, -nzv.col]
  
  
  
 
  d.pca <- prcomp(data.pca,
                  center = TRUE,
                  scale. = TRUE)
  
  # Create percent variance label for axis 
  percentVar <- d.pca$sdev^2 / sum( d.pca$sdev^2 )
  percentVar <- round(100 * percentVar)
  
  #Order groupings 
  plot.names <- factor(plot.names, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
  grouping.colours <- c("#40004b", "#762a83", "#c2a5cf", "#b8e186", "#4d9221", "#276419", "#fed976", "#fd8d3c", "#bd0026")
  
  (pcaplot <- ggplot(d.pca, aes(PC1, PC2, color = plot.names , fill = plot.names)) +
      coord_fixed() +
      geom_point(size=3) +
      scale_fill_manual(values = grouping.colours) +
      scale_colour_manual(values = grouping.colours) +
      #must manually apply the percentage variance in this case
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      #geom_polygon(data = hulls, alpha = 0.5) +
      stat_chull(geom = "polygon", alpha = 0.2) +
      #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
      #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
      stat_mean(size = 5) +
      theme_bw() +
      theme(aspect.ratio=1,
            legend.title = element_blank(),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title=element_text(size = 22),
            legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_rpoB_PC12.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm' , useDingbats=FALSE)

  ########################
  ######### PC23 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC2, PC3, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC2: ",percentVar[2],"% variance")) +
     ylab(paste0("PC3: ",percentVar[3],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_rpoB_PC23.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
  
  ########################
  ######### PC34 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC3, PC4, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC3: ",percentVar[3],"% variance")) +
     ylab(paste0("PC4: ",percentVar[4],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_rpoB_PC34.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
}  


#================================================#
#                    recA                        #
#================================================#
{
  #Order 
  reca$index <- as.numeric(row.names(reca))
  reca <- reca[order(reca$index), ]
  samples <- reca$index
  reca$index <- NULL
  
  #Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  
  # Create plot column for colouring
  plot.names <- d2$Grouping
  
  data.pca <- reca
  
  # Replace 0s with 0.00001
  data.pca[data.pca == 0] <- 0.0000000001
  
  # Log transform data for pca
  # Change end row value to be the ncolumns
  x = ncol(data.pca)
  data.pca[ ,1:x] <- log10(data.pca[ ,1:x])
  
  # Remove 0 variance using caret for when Aarhus plots are removed previously
  NZV <- nearZeroVar(data.pca, saveMetrics = TRUE)
  nzv.col <- nearZeroVar(data.pca)
  if(length(nzv.col) > 0) data.pca <- data.pca[, -nzv.col]
  
  
  
  
  d.pca <- prcomp(data.pca,
                  center = TRUE,
                  scale. = TRUE)
  
  # Create percent variance label for axis 
  percentVar <- d.pca$sdev^2 / sum( d.pca$sdev^2 )
  percentVar <- round(100 * percentVar)
  
  #Order groupings 
  plot.names <- factor(plot.names, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
  grouping.colours <- c("#40004b", "#762a83", "#c2a5cf", "#b8e186", "#4d9221", "#276419", "#fed976", "#fd8d3c", "#bd0026")  
  
  (pcaplot <- ggplot(d.pca, aes(PC1, PC2, color = plot.names , fill = plot.names)) +
      coord_fixed() +
      geom_point(size=3) +       
      scale_fill_manual(values = grouping.colours) +       
      scale_colour_manual(values = grouping.colours) +
      #must manually apply the percentage variance in this case
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      #geom_polygon(data = hulls, alpha = 0.5) +
      stat_chull(geom = "polygon", alpha = 0.2) +
      #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
      #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
      stat_mean(size = 5) +
      theme_bw() +
      theme(aspect.ratio=1,
            legend.title = element_blank(),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title=element_text(size = 22),
            legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_recA_PC12.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm' , useDingbats=FALSE)
  
  ########################
  ######### PC23 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC2, PC3, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC2: ",percentVar[2],"% variance")) +
     ylab(paste0("PC3: ",percentVar[3],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_recA_PC23.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
  
  ########################
  ######### PC34 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC3, PC4, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC3: ",percentVar[3],"% variance")) +
     ylab(paste0("PC4: ",percentVar[4],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_recA_PC34.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
}  


#================================================#
#                    nodA                        #
#================================================#
{
  #Order 
  noda$index <- as.numeric(row.names(noda))
  noda <- noda[order(noda$index), ]
  samples <- noda$index
  noda$index <- NULL
  
  #Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  
  # Create plot column for colouring
  plot.names <- d2$Grouping
  
  data.pca <- noda
  
  # Replace 0s with 0.00001
  data.pca[data.pca == 0] <- 0.0000000001
  
  # Log transform data for pca
  # Change end row value to be the ncolumns
  x = ncol(data.pca)
  data.pca[ ,1:x] <- log10(data.pca[ ,1:x])
  
  # Remove 0 variance using caret for when Aarhus plots are removed previously
  NZV <- nearZeroVar(data.pca, saveMetrics = TRUE)
  nzv.col <- nearZeroVar(data.pca)
  if(length(nzv.col) > 0) data.pca <- data.pca[, -nzv.col]
  
  
  
  
  d.pca <- prcomp(data.pca,
                  center = TRUE,
                  scale. = TRUE)
  
  # Create percent variance label for axis 
  percentVar <- d.pca$sdev^2 / sum( d.pca$sdev^2 )
  percentVar <- round(100 * percentVar)
  
  #Order groupings 
  plot.names <- factor(plot.names, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
  grouping.colours <- c("#40004b", "#762a83", "#c2a5cf", "#b8e186", "#4d9221", "#276419", "#fed976", "#fd8d3c", "#bd0026")  
  
  (pcaplot <- ggplot(d.pca, aes(PC1, PC2, color = plot.names , fill = plot.names)) +
      coord_fixed() +
      geom_point(size=3) +       
      scale_fill_manual(values = grouping.colours) +       
      scale_colour_manual(values = grouping.colours) +
      #must manually apply the percentage variance in this case
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      #geom_polygon(data = hulls, alpha = 0.5) +
      stat_chull(geom = "polygon", alpha = 0.2) +
      #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
      #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
      stat_mean(size = 5) +
      theme_bw() +
      theme(aspect.ratio=1,
            legend.title = element_blank(),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title=element_text(size = 22),
            legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_nodA_PC12.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm' , useDingbats=FALSE)
  
  ########################
  ######### PC23 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC2, PC3, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC2: ",percentVar[2],"% variance")) +
     ylab(paste0("PC3: ",percentVar[3],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_nodA_PC23.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
  
  ########################
  ######### PC34 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC3, PC4, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC3: ",percentVar[3],"% variance")) +
     ylab(paste0("PC4: ",percentVar[4],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_nodA_PC34.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
}


#================================================#
#                    nodD                        #
#================================================#
{
  #Order 
  nodd$index <- as.numeric(row.names(nodd))
  nodd <- nodd[order(nodd$index), ]
  samples <- nodd$index
  nodd$index <- NULL
  
  #Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  
  # Create plot column for colouring
  plot.names <- d2$Grouping
  
  data.pca <- nodd
  
  # Replace 0s with 0.00001
  data.pca[data.pca == 0] <- 0.0000000001
  
  # Log transform data for pca
  # Change end row value to be the ncolumns
  x = ncol(data.pca)
  data.pca[ ,1:x] <- log10(data.pca[ ,1:x])
  
  # Remove 0 variance using caret for when Aarhus plots are removed previously
  NZV <- nearZeroVar(data.pca, saveMetrics = TRUE)
  nzv.col <- nearZeroVar(data.pca)
  if(length(nzv.col) > 0) data.pca <- data.pca[, -nzv.col]
  
  
  
  
  d.pca <- prcomp(data.pca,
                  center = TRUE,
                  scale. = TRUE)
  
  # Create percent variance label for axis 
  percentVar <- d.pca$sdev^2 / sum( d.pca$sdev^2 )
  percentVar <- round(100 * percentVar)
  
  #Order groupings 
  plot.names <- factor(plot.names, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
  grouping.colours <- c("#40004b", "#762a83", "#c2a5cf", "#b8e186", "#4d9221", "#276419", "#fed976", "#fd8d3c", "#bd0026")
  
  (pcaplot <- ggplot(d.pca, aes(PC1, PC2, color = plot.names , fill = plot.names)) +
      coord_fixed() +
      geom_point(size=3) +       
      scale_fill_manual(values = grouping.colours) +       
      scale_colour_manual(values = grouping.colours) +
      #must manually apply the percentage variance in this case
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      #geom_polygon(data = hulls, alpha = 0.5) +
      stat_chull(geom = "polygon", alpha = 0.2) +
      #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
      #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
      stat_mean(size = 5) +
      theme_bw() +
      theme(aspect.ratio=1,
            legend.title = element_blank(),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title=element_text(size = 22),
            legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_nodD_PC12.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm' , useDingbats=FALSE)
  
  ########################
  ######### PC23 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC2, PC3, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC2: ",percentVar[2],"% variance")) +
     ylab(paste0("PC3: ",percentVar[3],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_nodD_PC23.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
  
  ########################
  ######### PC34 #########
  ########################
  (pcaplot <- ggplot(d.pca, aes(PC3, PC4, color = plot.names , fill = plot.names)) +
     coord_fixed() +
     geom_point(size=3) +       
     scale_fill_manual(values = grouping.colours) +       
     scale_colour_manual(values = grouping.colours) +
     #must manually apply the percentage variance in this case
     xlab(paste0("PC3: ",percentVar[3],"% variance")) +
     ylab(paste0("PC4: ",percentVar[4],"% variance")) + 
     #geom_polygon(data = hulls, alpha = 0.5) +
     stat_chull(geom = "polygon", alpha = 0.2) +
     #stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95) + #if you want to look at ellipses of t-distibuted data 95%. 
     #stat_conf_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, bary = TRUE) + #if you want to look at ellipses of mean 95%. bary = false is the same as stat_ellipse(type = "euclid", level - 0.95).
     stat_mean(size = 5) +
     theme_bw() +
     theme(aspect.ratio=1,
           legend.title = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           axis.title=element_text(size = 22),
           legend.text = element_text(size = 20)))
  
  
  ggsave('PCA_nodD_PC34.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')
  
}


## Colour scale with blue-yellow-red


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


