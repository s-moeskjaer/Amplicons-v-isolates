#install.packages("vegan")
#install.packages("MASS")
#install.packages("pheatmap")
library(vegan)
library(MASS)
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)
library(dplyr)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(multcompView)


setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/')

#==========================================#
#                                          #
#            Load all data                 #
#                                          #
#==========================================#
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



#==========================================#
#                                          #
#            Allele abundance              #
#                                          #
#==========================================#

########## rpoB ###########
#Order 
rpob$index <- as.numeric(row.names(rpob))
rpob <- rpob[order(rpob$index), ]
samples <- rpob$index
rpob$index <- NULL

#Add annotation row based on each dataset
d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
#reorder genospecies and associate a colour 
d2$Origin <- factor(d2$Origin, levels = c("DKO", "DK", "F", "UK"))
row.names(d2) <- samples

allele_n <- rowSums(rpob != 0)
d2$allele_n <- allele_n
ggplot(d2, aes(x=Sample, y=allele_n, color = Origin)) +
  geom_point() 


########## recA ###########
#Order 
reca$index <- as.numeric(row.names(reca))
reca <- reca[order(reca$index), ]
samples <- reca$index
reca$index <- NULL

#Add annotation row based on each dataset
d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
#reorder genospecies and associate a colour 
d2$Origin <- factor(d2$Origin, levels = c("DKO", "DK", "F", "UK"))
row.names(d2) <- samples

allele_n <- rowSums(reca != 0)
d2$allele_n <- allele_n
ggplot(d2, aes(x=Sample, y=allele_n, color = Origin)) +
  geom_point() 


########## nodA ###########
#Order 
noda$index <- as.numeric(row.names(noda))
noda <- noda[order(noda$index), ]
samples <- noda$index
noda$index <- NULL

#Add annotation row based on each dataset
d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
#reorder genospecies and associate a colour 
d2$Origin <- factor(d2$Origin, levels = c("DKO", "DK", "F", "UK"))
row.names(d2) <- samples

allele_n <- rowSums(noda != 0)
d2$allele_n <- allele_n
ggplot(d2, aes(x=Sample, y=allele_n, color = Origin)) +
  geom_point() 



########## nodD ###########
#Order 
nodd$index <- as.numeric(row.names(nodd))
nodd <- nodd[order(nodd$index), ]
samples <- nodd$index
nodd$index <- NULL

#Add annotation row based on each dataset
d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
#reorder genospecies and associate a colour 
d2$Origin <- factor(d2$Origin, levels = c("DKO", "DK", "F", "UK"))
row.names(d2) <- samples

allele_n <- rowSums(nodd != 0)
d2$allele_n <- allele_n
ggplot(d2, aes(x=Sample, y=allele_n, color = Origin)) +
  geom_point() 



#==========================================#
#                                          #
#                 Combined                 #
#                                          #
#==========================================#

#Order 
data_norm$index <- as.numeric(row.names(data_norm))
data_norm <- data_norm[order(data_norm$index), ]
samples <- data_norm$index
data_norm$index <- NULL
#Add annotation row based on each dataset
d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]

#reorder genospecies and associate a colour 
#d2$Origin <- factor(d2$Origin, levels = c("DKO", "DK", "F", "UK"))

#reorder groupings and associate a colour 
d2$Grouping <- factor(d2$Grouping, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))

row.names(d2) <- samples


d2$nodA <- rowSums(data_norm[,30:51] != 0)
d2$nodD <- rowSums(data_norm[,52:72] != 0)
d2$recA <- rowSums(data_norm[,17:29] != 0)
d2$rpoB <- rowSums(data_norm[,1:16] != 0)

d3 <- d2[,9:12]
d3$Grouping <- d2$Grouping
forplot <- melt(d3, id.vars="Grouping")

###### Group by origin
plot1 <- ggplot(forplot, aes(x=Grouping, y=value, fill = variable)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.size=0.5,color = 'black') +
  theme_bw() +
  labs(x="Origin", y = "n unique alleles") +
  #scale_fill_manual(values=c("#b2182b", "#f4a582", "#e0e0e0", "#4d4d4d")) # grey-red colours
  scale_fill_manual(values=c("#01665e", "#c7eae5", "#dfc27d", "#8c510a")) + # petroleum-brown colours
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size = 22),
        legend.text = element_text(size = 20))

plot1
ggsave('allele_number.pdf', plot = plot1, width = 25, height = 20, unit = 'cm')

####### Group by gene 
plot2 <- ggplot(forplot, aes(x=variable, y=value, fill = Grouping)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.size=0.5,color = 'black') +
  theme_bw() +
  labs(x="Origin", y = "n unique alleles") +
  #scale_fill_manual(values=c("#b2182b", "#f4a582", "#e0e0e0", "#4d4d4d")) # grey-red colours
  #scale_fill_manual(values=c("#0f1027", "#5c2e88", "#ba4f8f", "#f68f75")) + # petroleum-brown colours
  scale_fill_manual(values=c("#40004b", "#762a83", "#c2a5cf", "#b8e186", "#4d9221", "#276419", "#fed976", "#fd8d3c", "#bd0026")) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size = 22),
        legend.text = element_text(size = 20))

plot2
ggsave('allele_number_grouping1.pdf', plot = plot2, width = 30, height = 15, unit = 'cm')



#=============#
#     rpoB    #
#=============#
d3 <- data.frame(d2$rpoB,d2$Grouping)
forplot <- melt(d3, id.vars="d2.Grouping")

####### Group by gene 
plot2 <- ggplot(forplot, aes(x=d2.Grouping, y=value, fill = d2.Grouping)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.size=0.5,color = 'black') +
  theme_bw() +
  labs(x="Origin", y = "n unique alleles") +
  #scale_fill_manual(values=c("#b2182b", "#f4a582", "#e0e0e0", "#4d4d4d")) # grey-red colours
  #scale_fill_manual(values=c("#0f1027", "#5c2e88", "#ba4f8f", "#f68f75")) + # petroleum-brown colours
  scale_fill_manual(values=c("#6e016b", "#88419d", "#8c6bb1", "#8c96c6", "#9ebcda", "#bfd3e6", "#ccebc5", "#4eb3d3", "#08589e")) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size = 22),
        legend.text = element_text(size = 20))

plot2
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/allele_counts/rpoB_grouping.pdf', plot = plot2, width = 30, height = 15, unit = 'cm')

# ANOVA
fit <- aov(d2.rpoB ~ factor(d2.Grouping)  , data=d3)
results <- TukeyHSD(fit)
multcompLetters4(fit, results)



#=============#
#     recA    #
#=============#
d3 <- data.frame(d2$recA,d2$Grouping)
forplot <- melt(d3, id.vars="d2.Grouping")

####### Group by gene 
plot2 <- ggplot(forplot, aes(x=d2.Grouping, y=value, fill = d2.Grouping)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.size=0.5,color = 'black') +
  theme_bw() +
  labs(x="Origin", y = "n unique alleles") +
  #scale_fill_manual(values=c("#b2182b", "#f4a582", "#e0e0e0", "#4d4d4d")) # grey-red colours
  #scale_fill_manual(values=c("#0f1027", "#5c2e88", "#ba4f8f", "#f68f75")) + # petroleum-brown colours
  scale_fill_manual(values=c("#6e016b", "#88419d", "#8c6bb1", "#8c96c6", "#9ebcda", "#bfd3e6", "#ccebc5", "#4eb3d3", "#08589e")) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size = 22),
        legend.text = element_text(size = 20))

plot2
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/allele_counts/recA_grouping.pdf', plot = plot2, width = 30, height = 15, unit = 'cm')

# ANOVA
fit <- aov(d2.recA ~ factor(d2.Grouping)  , data=d3)
results <- TukeyHSD(fit)
multcompLetters4(fit, results)



#=============#
#     nodA    #
#=============#
d3 <- data.frame(d2$nodA,d2$Grouping)
forplot <- melt(d3, id.vars="d2.Grouping")

####### Group by gene 
plot2 <- ggplot(forplot, aes(x=d2.Grouping, y=value, fill = d2.Grouping)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.size=0.5,color = 'black') +
  theme_bw() +
  labs(x="Origin", y = "n unique alleles") +
  #scale_fill_manual(values=c("#b2182b", "#f4a582", "#e0e0e0", "#4d4d4d")) # grey-red colours
  #scale_fill_manual(values=c("#0f1027", "#5c2e88", "#ba4f8f", "#f68f75")) + # petroleum-brown colours
  scale_fill_manual(values=c("#6e016b", "#88419d", "#8c6bb1", "#8c96c6", "#9ebcda", "#bfd3e6", "#ccebc5", "#4eb3d3", "#08589e")) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size = 22),
        legend.text = element_text(size = 20))

plot2
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/allele_counts/nodA_grouping.pdf', plot = plot2, width = 30, height = 15, unit = 'cm')

# ANOVA
fit <- aov(d2.nodA ~ factor(d2.Grouping)  , data=d3)
results <- TukeyHSD(fit)
multcompLetters4(fit, results)


#=============#
#     nodD    #
#=============#
d3 <- data.frame(d2$nodD,d2$Grouping)
forplot <- melt(d3, id.vars="d2.Grouping")

####### Group by gene 
plot2 <- ggplot(forplot, aes(x=d2.Grouping, y=value, fill = d2.Grouping)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.size=0.5,color = 'black') +
  theme_bw() +
  labs(x="Origin", y = "n unique alleles") +
  #scale_fill_manual(values=c("#b2182b", "#f4a582", "#e0e0e0", "#4d4d4d")) # grey-red colours
  #scale_fill_manual(values=c("#0f1027", "#5c2e88", "#ba4f8f", "#f68f75")) + # petroleum-brown colours
  scale_fill_manual(values=c("#6e016b", "#88419d", "#8c6bb1", "#8c96c6", "#9ebcda", "#bfd3e6", "#ccebc5", "#4eb3d3", "#08589e")) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size = 22),
        legend.text = element_text(size = 20))

plot2
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/allele_counts/nodD_grouping.pdf', plot = plot2, width = 30, height = 15, unit = 'cm')

# ANOVA
fit <- aov(d2.nodD ~ factor(d2.Grouping)  , data=d3)
results <- TukeyHSD(fit)
multcompLetters4(fit, results)
