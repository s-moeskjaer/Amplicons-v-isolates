# https://popgen.nescent.org/DifferentiationSNP.html
# Verbose code explanation at the bottom of document
library("adegenet")
library("hierfstat")
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/')

# meta data 
envir <- read.table('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/zz_sample_overview_latlong.csv', header = T, sep = ";")
names(envir)[names(envir)=="X...order"] <- "order"

data_soil = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Table_S1.csv', sep = ';')

#==================================================================================================================#
#                                                                                                                  #
#                                                    GROUPINGS                                                     #
#                                                                                                                  #
#==================================================================================================================#
#=============================#
#             rpoB            #
#=============================#
df <- read.table("rpoB_input_norm.txt", header = TRUE, check.names = FALSE)   


# ========================================== #


ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$grouping) # use later with adegenet (population labels)
county <- df$sample 

#To convert df to a genind object (adegenet), the input should only contain genotypes
x = ncol(df)
locus   <- df[, 5:x] 
df1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

# Pairwise Fst
#genet.dist(df1, method = "WC84")
fst_matrix <- as.matrix(genet.dist(df1, method = "WC84"))
fst_matrix1 <- fst_matrix[c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"),]
fst_matrix2 <- fst_matrix1[,c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK")]

#Add annotation row based on each dataset
samples <- row.names(fst_matrix2)

#reorder genospecies and associate a colour 
d2 <- data.frame(row.names(fst_matrix2), samples)
row.names(d2) <- row.names(fst_matrix2)
d2$samples <- factor(d2$samples, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
d3 <- d2[2]

my_palette = rev(brewer.pal(9, "Greys"))
#Colour list must be named after column in d2
anno_colors3 <- list(Grouping = c(DKO_1 = "#40004b", DKO_2 = "#762a83", DKO_3 = "#c2a5cf", DKO_4 = "#b8e186", DKO_5 = "#4d9221", DKO_6 = "#276419", DK = "#fed976", F = "#fd8d3c", UK = "#bd0026"))
#Plot
pheatmap(fst_matrix2,
         cluster_rows = F, cluster_cols = F,
         #gaps_row = c(length(which(d3$Origin=="DKO")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")) + length(which(d3$Origin=="F"))),
         #gaps_col = c(length(which(d3$Origin=="DKO")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")) + length(which(d3$Origin=="F"))),
         annotation_row = d3,
         annotation_col = d3,
         annotation_colors = anno_colors3,
         show_colnames     = FALSE,
         show_rownames     = FALSE,
         border_color = NA,
         col = (my_palette),
         filename = "/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/hierf_rpoB_fst_group.pdf")


#=============================#
#             recA            #
#=============================#
df <- read.table("recA_input_norm.txt", header = TRUE, check.names = FALSE)   


ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$grouping) # use later with adegenet (population labels)
county <- df$sample 

#To convert df to a ???genind??? object (adegenet), the input should only contain genotypes
x = ncol(df)
locus   <- df[, 5:x] 
df1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

# Pairwise Fst
#genet.dist(df1, method = "WC84")
fst_matrix <- as.matrix(genet.dist(df1, method = "WC84"))

#Add annotation row based on each dataset
samples <- row.names(fst_sample_matrix)
d2 <- envir[which(envir$Sample %in% samples),]
#d2 <- d2[order(d2$order), ]
#reorder genospecies and associate a colour 
d2$Grouping <- factor(d2$Grouping, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
row.names(d2) <- row.names(fst_sample_matrix)
d3 <- d2[4]
d2 <- d2[8]

my_palette = rev(brewer.pal(9, "Greys"))
#Colour list must be named after column in d2
anno_colors3 <- list(Grouping = c(DKO_1 = "#40004b", DKO_2 = "#762a83", DKO_3 = "#c2a5cf", DKO_4 = "#b8e186", DKO_5 = "#4d9221", DKO_6 = "#276419", DK = "#fed976", F = "#fd8d3c", UK = "#bd0026"))
#Plot
pheatmap(fst_sample_matrix,
         #cluster_rows = F, cluster_cols = F,
         #gaps_row = c(length(which(d3$Origin=="DKO")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")) + length(which(d3$Origin=="F"))),
         #gaps_col = c(length(which(d3$Origin=="DKO")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")) + length(which(d3$Origin=="F"))),
         annotation_row = d2,
         annotation_col = d2,
         annotation_colors = anno_colors3,
         show_colnames     = FALSE,
         show_rownames     = FALSE,
         border_color = NA,
         col = (my_palette),
         filename = "/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/hierf_rpoB_fst_all.pdf")

#=============================#
#             nodA            #
#=============================#
df <- read.table("nodA_input_norm.txt", header = TRUE, check.names = FALSE)   
df <- df[,-7]

#DKO only 
y = nrow(df)
df <- df[-3414:-y,]

ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$grouping) # use later with adegenet (population labels)
county <- df$sample 

#To convert df to a ???genind??? object (adegenet), the input should only contain genotypes
x = ncol(df)
locus   <- df[, 5:x] 
df1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

# Pairwise Fst
#genet.dist(df1, method = "WC84")
fst_matrix <- as.matrix(genet.dist(df1, method = "WC84"))

# GPS coordinates
gps_groups <- read.csv("grouping_gps.csv", header = T, sep = ";")
samples <- row.names(fst_matrix)
d2 <- gps_groups[which(gps_groups$Grouping %in% samples),]
## Sampling site points 
gps <- data.frame(d2$Long,d2$Lat) 
distmat <- distm(gps, gps, fun = distVincentyEllipsoid)

# plotting using ggplot 
diag(distmat) = NA 
diag(fst_matrix) = NA
x <- as.vector( t( distmat ))
y <- as.vector(t(fst_matrix))
correlation(x,y) 

df = data.frame(x,y)
p <- ggplot( df, aes( x = x/1000, y = y)) + 
  geom_point(colour = "#7a0177") +
  geom_smooth(colour = "#f768a1", fill = "#fcc5c0", method = 'lm') +
  scale_y_continuous(limits = c(-0.02, 0.13)) +
  #ggtitle( "nodA, cor = 0.6316436, p-value = 0.0001814762") +
  labs( x = "Pairwise distance (km)", y = "Fst" ) + 
  theme_bw() 

p
#ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_group_noda.pdf', plot = p, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_group_nodatest.pdf', plot = p, width = 7, height = 7, unit = 'cm', useDingbats = FALSE)


#=============================#
#             nodD            #
#=============================#
df <- read.table("nodD_input_norm.txt", header = TRUE, check.names = FALSE)   

#DKO only 
y = nrow(df)
df <- df[-3700:-y,]

ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$grouping) # use later with adegenet (population labels)
county <- df$sample 

#To convert df to a ???genind??? object (adegenet), the input should only contain genotypes
x = ncol(df)
locus   <- df[, 5:x] 
df1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

# Pairwise Fst
#genet.dist(df1, method = "WC84")
fst_matrix <- as.matrix(genet.dist(df1, method = "WC84"))

# GPS coordinates
gps_groups <- read.csv("grouping_gps.csv", header = T, sep = ";")
samples <- row.names(fst_matrix)
d2 <- gps_groups[which(gps_groups$Grouping %in% samples),]
## Sampling site points 
gps <- data.frame(d2$Long,d2$Lat) 
distmat <- distm(gps, gps, fun = distVincentyEllipsoid)

# plotting using ggplot 
diag(distmat) = NA 
diag(fst_matrix) = NA
x <- as.vector( t( distmat ))
y <- as.vector(t(fst_matrix))
correlation(x,y) 

df = data.frame(x,y)
p <- ggplot( df, aes( x = x/1000, y = y)) + 
  geom_point(colour = "#7a0177") +
  geom_smooth(colour = "#f768a1", fill = "#fcc5c0", method = 'lm') +
  scale_y_continuous(limits = c(-0.02, 0.13)) +
  #ggtitle( "nodD, cor = 0.4204567, p-value = 0.02069601") +
  labs( x = "Pairwise distance (km)", y = "Fst" ) + 
  theme_bw() 

p
#ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_group_nodd.pdf', plot = p, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_group_noddtest.pdf', plot = p, width = 7, height = 7, unit = 'cm', useDingbats = FALSE)





#==================================================================================================================#
#                                                                                                                  #
#                                                  PER SAMPLE                                                      #
#                                                                                                                  #
#==================================================================================================================#

#=============================#
#             rpoB            #
#=============================#
df <- read.table("rpoB_input_norm.txt", header = TRUE, check.names = FALSE)   

ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$sample) # use later with adegenet (population labels)

#To convert df to a ???genind??? object (adegenet), the input should only contain genotypes
x = ncol(df)
locus   <- df[, 5:x] 
df1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

# Pairwise Fst
#genet.dist(df1, method = "WC84")
fst_sample_matrix <- as.matrix(genet.dist(df1, method = "WC84"))
fst_sample_matrix[fst_sample_matrix < 0] <- 0

#Add annotation row based on each dataset
samples <- row.names(fst_sample_matrix)
d2 <- envir[which(envir$Sample %in% samples),]
#d2 <- d2[order(d2$order), ]
#reorder genospecies and associate a colour 
d2$Grouping <- factor(d2$Grouping, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
row.names(d2) <- row.names(fst_sample_matrix)
d3 <- d2[4]
d2 <- d2[8]

my_palette = rev(brewer.pal(9, "Greys"))
#Colour list must be named after column in d2
anno_colors3 <- list(Grouping = c(DKO_1 = "#40004b", DKO_2 = "#762a83", DKO_3 = "#c2a5cf", DKO_4 = "#b8e186", DKO_5 = "#4d9221", DKO_6 = "#276419", DK = "#fed976", F = "#fd8d3c", UK = "#bd0026"))
#Plot
pheatmap(fst_sample_matrix,
         #cluster_rows = F, cluster_cols = F,
         #gaps_row = c(length(which(d3$Origin=="DKO")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")) + length(which(d3$Origin=="F"))),
         #gaps_col = c(length(which(d3$Origin=="DKO")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")) + length(which(d3$Origin=="F"))),
         annotation_row = d2,
         annotation_col = d2,
         annotation_colors = anno_colors3,
         show_colnames     = FALSE,
         show_rownames     = FALSE,
         border_color = NA,
         col = (my_palette),
         filename = "/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/hierf_rpoB_fst_all.pdf")

#=============================#
#             recA            #
#=============================#
df <- read.table("recA_input_norm.txt", header = TRUE, check.names = FALSE)   

ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$sample) # use later with adegenet (population labels)

#To convert df to a ???genind??? object (adegenet), the input should only contain genotypes
x = ncol(df)
locus   <- df[, 5:x] 
df1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

# Pairwise Fst
#genet.dist(df1, method = "WC84")
fst_sample_matrix <- as.matrix(genet.dist(df1, method = "WC84"))
fst_sample_matrix[fst_sample_matrix < 0] <- 0

#Add annotation row based on each dataset
samples <- row.names(fst_sample_matrix)
d2 <- envir[which(envir$Sample %in% samples),]
#d2 <- d2[order(d2$order), ]
#reorder genospecies and associate a colour 
d2$Grouping <- factor(d2$Grouping, levels = c("DKO_1", "DKO_2", "DKO_3", "DKO_4", "DKO_5", "DKO_6", "DK", "F", "UK"))
row.names(d2) <- row.names(fst_sample_matrix)
d3 <- d2[4]
d2 <- d2[8]

my_palette = rev(brewer.pal(9, "Greys"))
#Colour list must be named after column in d2
anno_colors3 <- list(Grouping = c(DKO_1 = "#40004b", DKO_2 = "#762a83", DKO_3 = "#c2a5cf", DKO_4 = "#b8e186", DKO_5 = "#4d9221", DKO_6 = "#276419", DK = "#fed976", F = "#fd8d3c", UK = "#bd0026"))
#Plot
pheatmap(fst_sample_matrix,
         #cluster_rows = F,
         cluster_rows = F, cluster_cols = F,
         #gaps_row = c(length(which(d3$Origin=="DKO")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")) + length(which(d3$Origin=="F"))),
         #gaps_col = c(length(which(d3$Origin=="DKO")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")), length(which(d3$Origin=="DKO")) + length(which(d3$Origin=="DK")) + length(which(d3$Origin=="F"))),
         annotation_row = d2,
         annotation_col = d2,
         annotation_colors = anno_colors3,
         show_colnames     = FALSE,
         show_rownames     = FALSE,
         border_color = NA,
         col = (my_palette),
         filename = "/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/hierf_recA_fst_all.pdf")


#=============================#
#             nodA            #
#=============================#
df <- read.table("nodA_input_norm.txt", header = TRUE, check.names = FALSE)   
df <- df[,-7]

#DKO only 
y = nrow(df)
df <- df[-3414:-y,]

ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$sample) # use later with adegenet (population labels)

#To convert df to a ???genind??? object (adegenet), the input should only contain genotypes
x = ncol(df)
locus   <- df[, 5:x] 
df1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

# Pairwise Fst
#genet.dist(df1, method = "WC84")
fst_sample_matrix <- as.matrix(genet.dist(df1, method = "WC84"))
fst_sample_matrix[fst_sample_matrix < 0] <- 0

# GPS coordinates
samples <- row.names(fst_sample_matrix)
d2 <- envir[which(envir$Sample %in% samples),]
## Sampling site points 
gps <- data.frame(d2$Long,d2$Lat) 
distmat <- distm(gps, gps, fun = distVincentyEllipsoid)

# plotting using ggplot 
diag(distmat) = NA 
diag(fst_sample_matrix) = NA
x <- as.vector( t( distmat ))
y <- as.vector(t(fst_sample_matrix))
correlation(x,y)  

df = data.frame(x,y)
p <- ggplot( df, aes( x = x/1000, y = y)) + 
  geom_point(colour = "#7a0177") +
  geom_smooth(colour = "#f768a1", fill = "#fcc5c0", method = 'lm') +
  scale_y_continuous(limits = c(0, 0.22)) +
  #ggtitle( "nodA, cor = 0.02033695, p-value = 0.4961734") +
  labs( x = "Pairwise distance (km)", y = "Fst" ) + 
  theme_bw() 
p
#ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_nodA.pdf', plot = p, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_nodAtest.pdf', plot = p, width = 7, height = 7, unit = 'cm', useDingbats = FALSE)

#=============================#
#             nodD            #
#=============================#
df <- read.table("nodD_input_norm.txt", header = TRUE, check.names = FALSE)   

#DKO only 
y = nrow(df)
df <- df[-3700:-y,]

ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$sample) # use later with adegenet (population labels)

#To convert df to a ???genind??? object (adegenet), the input should only contain genotypes
x = ncol(df)
locus   <- df[, 5:x] 
df1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

# Pairwise Fst
#genet.dist(df1, method = "WC84")
fst_sample_matrix <- as.matrix(genet.dist(df1, method = "WC84"))
fst_sample_matrix[fst_sample_matrix < 0] <- 0

# GPS coordinates
samples <- row.names(fst_sample_matrix)
d2 <- envir[which(envir$Sample %in% samples),]
## Sampling site points 
gps <- data.frame(d2$Long,d2$Lat) 
distmat <- distm(gps, gps, fun = distVincentyEllipsoid)

# plotting using ggplot 
diag(distmat) = NA 
diag(fst_sample_matrix) = NA
x <- as.vector( t( distmat ))
y <- as.vector(t(fst_sample_matrix))
correlation(x,y)  

df = data.frame(x,y)
p <- ggplot( df, aes( x = x/1000, y = y)) + 
  geom_point(colour = "#7a0177") +
  geom_smooth(colour = "#f768a1", fill = "#fcc5c0", method = 'lm') +
  scale_y_continuous(limits = c(0, 0.22)) +
  #ggtitle( "nodD, cor = -0.007324366, p-value = 0.7894177") +
  labs( x = "Pairwise distance (km)", y = "Fst" ) + 
  theme_bw() 
p
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_nodD.pdf', plot = p, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_nodDtest.pdf', plot = p, width = 7, height = 7, unit = 'cm', useDingbats = FALSE)





#==================================================================================================================#
#                                                                                                                  #
#                                              CENTROID GROUPINGS                                                  #
#                                                                                                                  #
#==================================================================================================================#
d2 <- subset(envir, envir$Grouping == "DKO_6")
gps <- data.frame(d2$Long,d2$Lat) 
centroid(gps)

#==================================================================================================================#
#                                                                                                                  #
#                                                     VERBOSE                                                      #
#                                                                                                                  #
#==================================================================================================================#

df <- read.table("nodA_DKO_Fst.txt", header = TRUE, check.names = FALSE)   
dim(df) 

ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$group) # use later with adegenet (population labels)
county <- df$sample 

#To convert df to a ???genind??? object (adegenet), the input should only contain genotypes
locus   <- df[, 5:29] 
df1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = "")
df1

df2 <- genind2hierfstat(df1) 

basic.stats(df1) # Fst following Nei (1987) on genind object

wc(df1) # Weir and Cockerham's estimate


#=======================================#
#         Hierarchical Fst tests        #
#        =AMOVA for SNP dataset         #
#=======================================#
#The function varcomp.glob() produces a Hierarchical Fst (=AMOVA for SNPs or bi-allelic markers) 
#It is possible to make permutations on the different levels: 
#The function test.g() tests the effect of the population on genetic differentiation. 
#Individuals are randomly permuted among states. 
#The states influence genetic differentiation at a 5% level. 
#With the function test.between(), the counties are permuted among states. 
#The states influence significantly genetic structuring.
loci <- df2[, -1] # Remove the population column
varcomp.glob(levels = data.frame(population, county), loci, diploid = TRUE) 

test.g(loci, level = population) 

test.between(loci, test.lev = population, rand.unit = county, nperm = 100) 


#=============================#
#         Pairwise Fst        #
#=============================#
genet.dist(df1, method = "WC84")

fst_matrix <- as.matrix(genet.dist(df1, method = "WC84"))

#================================#
#    Unsupervised clustering     #
#================================#
# using Kmeans and DAPC in adegenet 
set.seed(20160308) # Setting a seed for a consistent result
grp <- find.clusters(df1, max.n.clust = 10, n.pca = 25, choose.n.clust = FALSE) 
names(grp)

dapc1 <- dapc(df1, grp$grp, n.pca = 20, n.da = 6) 
scatter(dapc1) # plot of the group


