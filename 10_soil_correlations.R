#http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
library(reshape2)
library(plyr)
library(data.table)
library(corrgram)
library(corrplot)
library(ggplot2)
library(agricolae)
library(GGally)
library(RColorBrewer)
library(caret)


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
#                   rpoB                   #
#                                          #
#==========================================#
{
  # Add soil data for correlations 
  data_soil = read.csv('Table_S1.csv', sep = ';')
  names(data_soil)[names(data_soil)=="Field.of.origin"] <- "index"
  
  rpob$index <- as.numeric(row.names(rpob))
  rpob <- rpob[order(rpob$index), ]
  samples <- rpob$index
  
  
  # Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  rpob$index <- d2$ID
  
  
  
  #==========================================#
  #                    ALL                   #
  #==========================================#
  d1 <- merge(rpob,data_soil, by="index",all.x=FALSE)
  d2 <- data.frame(d1[2:17], d1[20:21], d1[24:33])
  plot <- corrgram(d2, order=TRUE, lower.panel=panel.shade,
           upper.panel=panel.pts, text.panel=panel.txt)
  
  # Correlation matrix
  M<-cor(d2)
  pdf("corplot_rpoB_all_unordered.pdf", bg = "white")
  corrplot(M, type="upper", tl.col="black")
  dev.off()
  pdf("corplot_rpoB_all_ordered.pdf", bg = "white")
  corrplot(M, type="upper", order="hclust", tl.col="black")
  dev.off()
  
  #==========================================#
  #                    DKO                   #
  #==========================================#
  rpobDKO <- rpob[1:39,]
  d1 <- merge(rpobDKO,data_soil, by="index",all.x=FALSE)
  
  d2 <- data.frame(d1[2:17], d1[20:21], d1[24:33])
  d2$rpoBseq6 <- NULL
  d2$rpoBseq10 <- NULL
  # Add clay and silt to one column - slightly improves correlation to rpoBseq2, but not by much. Keep separate. 
  #d2$Clay_silt <- d2$Clay+d2$Silt
  #d2$Clay <- NULL
  #d2$Silt <- NULL
  
  plot <- corrgram(d2, order=F, lower.panel=panel.shade,
                   upper.panel=NULL, text.panel=panel.txt)
  
  # Correlations and scatter plots 
  correlation(d2$rpoBseq2,d2$Clay) # cor= 0.6140869, p-value = 5.025388e-05 
  correlation(d2$rpoBseq2,d2$Silt) # cor= 0.5876543, p-value = 0.0001048105 
  correlation(d2$rpoBseq2,d2$Coarse.sand) # cor= -0.5236968, p-value = 0.0007402673 
  correlation(d2$rpoBseq2,d2$Clay_silt) # cor= 0.6086513, p-value = 4.119025e-05 
  
  
  p1 <- ggplot(d2, aes(rpoBseq2, Clay)) +
    geom_point() +
    ggtitle( "cor =0.6140869, p-value = 4.119025e-05") 
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/soil_correlations/scatter_rpobseq2clay.pdf', plot = p1, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
  
  p2 <- ggplot(d2, aes(rpoBseq2, Silt)) +
    geom_point()+
    ggtitle( "cor =0.5876543, p-value = 0.0001048105") 
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/soil_correlations/scatter_rpobseq2silt.pdf', plot = p2, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
  
  p3 <- ggplot(d2, aes(rpoBseq2, Coarse.sand)) +
    geom_point()+
    ggtitle( "cor =-0.5236968, p-value = 0.0007402673") 
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/soil_correlations/scatter_rpobseq2coarsesand.pdf', plot = p3, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
  
  # Correlation matrix
  M<-cor(d2)
  pdf("corplot_rpoB_DKO_unordered.pdf", bg = "white")
  corrplot(M, type="upper", tl.col="black")
  dev.off()
  pdf("corplot_rpoB_DKO_ordered.pdf", bg = "white")
  corrplot(M, type="upper", order="hclust", tl.col="black")
  dev.off()
}

#==========================================#
#                                          #
#                   recA                   #
#                                          #
#==========================================#
{
  # Add soil data for correlations 
  data_soil = read.csv('Table_S1.csv', sep = ';')
  names(data_soil)[names(data_soil)=="Field.of.origin"] <- "index"
  
  reca$index <- as.numeric(row.names(reca))
  reca <- reca[order(reca$index), ]
  samples <- reca$index
  
  
  # Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  reca$index <- d2$ID
  
  
  
  #==========================================#
  #                    ALL                   #
  #==========================================#
  d1 <- merge(reca,data_soil, by="index",all.x=FALSE)
  d2 <- data.frame(d1[2:14], d1[17:18], d1[21:30])
  plot <- corrgram(d2, order=TRUE, lower.panel=panel.shade,
                   upper.panel=panel.pts, text.panel=panel.txt)
  
  # Correlation matrix
  M<-cor(d2)
  pdf("corplot_recA_all_unordered.pdf", bg = "white")
  corrplot(M, type="upper", tl.col="black")
  dev.off()
  pdf("corplot_recA_all_ordered.pdf", bg = "white")
  corrplot(M, type="upper", order="hclust", tl.col="black")
  dev.off()
  
  #==========================================#
  #                    DKO                   #
  #==========================================#
  recaDKO <- reca[1:46,]
  d1 <- merge(recaDKO,data_soil, by="index",all.x=FALSE)
  
  d2 <- data.frame(d1[2:14], d1[17:18], d1[21:30])
  d2$recAseq10 <- NULL
  
  plot <- corrgram(d2, order=F, lower.panel=panel.shade,
                   upper.panel=NULL, text.panel=panel.txt)
  
  # Correlations and scatter plots 
  correlation(d2$recAseq4,d2$Clay) # cor= 0.5106698, p-value = 0.0003380782 
  correlation(d2$recAseq6,d2$Silt) # cor= 0.1930562, p-value = 0.2038676 
  correlation(d2$recAseq1,d2$Clay) # cor= 0.1964362, p-value = 0.1959088 
  
  p1 <- ggplot(d2, aes(recAseq4, Clay)) +
    geom_point() +
    ggtitle( "cor =0.5106698, p-value = 0.0003380782") 
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/soil_correlations/scatter_recAseq4clay.pdf', plot = p1, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
  
  p2 <- ggplot(d2, aes(recAseq6, Silt)) +
    geom_point()+
    ggtitle( "cor =0.1930562, p-value = 0.2038676") 
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/soil_correlations/scatter_recAseq6silt.pdf', plot = p2, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
  
  p3 <- ggplot(d2, aes(recAseq1, Clay)) +
    geom_point()+
    ggtitle( "cor =0.1964362, p-value = 0.1959088") 
  ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/soil_correlations/scatter_recAseq1clay.pdf', plot = p3, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
  
  
  # Correlation matrix
  M<-cor(d2)
  pdf("corplot_recA_DKO_unordered.pdf", bg = "white")
  corrplot(M, type="upper", tl.col="black")
  dev.off()
  pdf("corplot_recA_DKO_ordered.pdf", bg = "white")
  corrplot(M, type="upper", order="hclust", tl.col="black")
  dev.off()
}



#==========================================#
#                                          #
#                   nodA                   #
#                                          #
#==========================================#
{
  # Add soil data for correlations 
  data_soil = read.csv('Table_S1.csv', sep = ';')
  names(data_soil)[names(data_soil)=="Field.of.origin"] <- "index"
  
  noda$index <- as.numeric(row.names(noda))
  noda <- noda[order(noda$index), ]
  samples <- noda$index
  
  
  # Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  noda$index <- d2$ID
  
  
  
  #==========================================#
  #                    ALL                   #
  #==========================================#
  d1 <- merge(noda,data_soil, by="index",all.x=FALSE)
  d2 <- data.frame(d1[2:22], d1[25:26], d1[29:38])
  plot <- corrgram(d2, order=TRUE, lower.panel=panel.shade,
                   upper.panel=panel.pts, text.panel=panel.txt)
  
  # Correlation matrix
  M<-cor(d2)
  pdf("corplot_nodA_all_unordered.pdf", bg = "white")
  corrplot(M, type="upper", tl.col="black")
  dev.off()
  pdf("corplot_nodA_all_ordered.pdf", bg = "white")
  corrplot(M, type="upper", order="hclust", tl.col="black")
  dev.off()
  
  #==========================================#
  #                    DKO                   #
  #==========================================#
  nodaDKO <- noda[1:34,]
  d1 <- merge(nodaDKO,data_soil, by="index",all.x=FALSE)
  
  d2 <- data.frame(d1[2:22], d1[25:26], d1[29:38])
  d2$nodAseq2 <- NULL
  d2$nodAseq11 <- NULL
  d2$nodAseq12 <- NULL
  d2$nodAseq19 <- NULL
  
  
  plot <- corrgram(d2, order=F, lower.panel=panel.shade,
                   upper.panel=NULL, text.panel=panel.txt)
  
  # Correlation matrix
  M<-cor(d2)
  pdf("corplot_nodA_DKO_unordered.pdf", bg = "white")
  corrplot(M, type="upper", tl.col="black")
  dev.off()
  pdf("corplot_nodA_DKO_ordered.pdf", bg = "white")
  corrplot(M, type="upper", order="hclust", tl.col="black")
  dev.off()
}


#==========================================#
#                                          #
#                   nodD                   #
#                                          #
#==========================================#
{
  # Add soil data for correlations 
  data_soil = read.csv('Table_S1.csv', sep = ';')
  names(data_soil)[names(data_soil)=="Field.of.origin"] <- "index"
  
  nodd$index <- as.numeric(row.names(nodd))
  nodd <- nodd[order(nodd$index), ]
  samples <- nodd$index
  
  
  # Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  nodd$index <- d2$ID
  
  
  
  #==========================================#
  #                    ALL                   #
  #==========================================#
  d1 <- merge(nodd,data_soil, by="index",all.x=FALSE)
  d2 <- data.frame(d1[2:18], d1[21:22], d1[25:34])
  plot <- corrgram(d2, order=TRUE, lower.panel=panel.shade,
                   upper.panel=panel.pts, text.panel=panel.txt)
  
  # Correlation matrix
  M<-cor(d2)
  pdf("corplot_nodD_all_unordered.pdf", bg = "white")
  corrplot(M, type="upper", tl.col="black")
  dev.off()
  pdf("corplot_nodD_all_ordered.pdf", bg = "white")
  corrplot(M, type="upper", order="hclust", tl.col="black")
  dev.off()
  
  #==========================================#
  #                    DKO                   #
  #==========================================#
  noddDKO <- nodd[1:37,]
  d1 <- merge(noddDKO,data_soil, by="index",all.x=FALSE)
  
  d2 <- data.frame(d1[2:18], d1[21:22], d1[25:34])
  d2$nodDseq1 <- NULL
  d2$nodDseq12 <- NULL
  d2$nodDseq14 <- NULL
  d2$nodDseq20 <- NULL
  
  plot <- corrgram(d2, order=F, lower.panel=panel.shade,
                   upper.panel=NULL, text.panel=panel.txt)
  
  # Correlation matrix
  M<-cor(d2)
  pdf("corplot_nodD_DKO_unordered.pdf", bg = "white")
  corrplot(M, type="upper", tl.col="black")
  dev.off()
  pdf("corplot_nodD_DKO_ordered.pdf", bg = "white")
  corrplot(M, type="upper", order="hclust", tl.col="black")
  dev.off()
}


#==========================================#
#                                          #
#                   ALL                    #
#                                          #
#==========================================#
{
  # Add soil data for correlations 
  data_soil = read.csv('Table_S1.csv', sep = ';')
  names(data_soil)[names(data_soil)=="Field.of.origin"] <- "index"
  
  samples <- row.names(data_norm)
  
  # Add annotation row based on each dataset
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  data_norm$index <- d2$ID
  
  
  
  #==========================================#
  #                    ALL                   #
  #==========================================#
  d1 <- merge(data_norm,data_soil, by="index",all.x=FALSE)
  d2 <- data.frame(d1[2:68], d1[71:73], d1[75:84])
  labs=colnames(d2)
  pdf("corplot_AllAll_ordered.pdf", bg = "white")
  corrgram(d2, order=TRUE, lower.panel=NULL,
                   upper.panel=panel.fill, text.panel=NULL,
                   outer.labels=list(top=list(labels=labs,cex=0.8),
                                     right=list(labels=labs,cex=0.8)))
           #col.regions=colorRampPalette(c("#5e4fa2", "#3288bd","#66c2a5", "#abdda4","#e6f598", "#fee08b","#fdae61", "#f46d43","#d53e4f","#9e0142")))
  dev.off()
  
  #==========================================#
  #                    DKO                   #
  #==========================================#
  data_norm2 <- data_norm[1:46,]
  d1 <- merge(data_norm2,data_soil, by="index",all.x=FALSE)
  d2 <- data.frame(d1[2:68], d1[71:73], d1[75:84])
  correlation(d2$rpoBseq2, d2$nodDseq2)
  plot(d2$rpoBseq2, d2$nodDseq2)

  # Remove 0 variance using caret for when Aarhus plots are removed previously
  NZV <- nearZeroVar(d2, saveMetrics = TRUE)
  nzv.col <- nearZeroVar(d2)
  if(length(nzv.col) > 0) d2 <- d2[, -nzv.col]
  
  labs=colnames(d2)
  pdf("corplot_AllAllDKO_ordered.pdf", bg = "white")
  corrgram(d2, order=TRUE, lower.panel=NULL,
           upper.panel=panel.fill, text.panel=NULL,
           outer.labels=list(top=list(labels=labs,cex=0.8),
                             right=list(labels=labs,cex=0.8)))
  dev.off()
}


