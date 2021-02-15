# https://popgen.nescent.org/DifferentiationSNP.html
# Verbose code explanation at the bottom of document
library("adegenet")
library("hierfstat")
library(geosphere)
library(agricolae)
library(ggplot2)
library(ade4)

setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/')

# meta data 
envir <- read.table('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/zz_sample_overview_latlong.csv', header = T, sep = ";")
names(envir)[names(envir)=="X...order"] <- "order"

#==================================================================================================================#
#                                                                                                                  #
#                                                    GROUPINGS                                                     #
#                                                                                                                  #
#==================================================================================================================#
#=============================#
#             rpoB            #
#=============================#
df <- read.table("rpoB_input_norm.txt", header = TRUE, check.names = FALSE)   

#DKO only 
y = nrow(df)
df <- df[-3914:-y,]

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

# ======== MANTEL TEST ======# 
geodist <- as.dist(distmat)
gendist <- as.dist(fst_matrix)
mantel.rtest(geodist, gendist, nrepet = 5000) # Observation: 0.7365868, Simulated p-value: 0.0019996  
#slope 
summary(lm(gendist~geodist))


# ======== Plotting ======# 
# plotting using ggplot 
diag(distmat) = NA 
diag(fst_matrix) = NA
x <- as.vector( t( distmat ))
y <- as.vector(t(fst_matrix))
correlation(x,y) 
xx <- x/1000
mlm = lm(y ~ xx) # slope = 0.000631
str(summary(mlm)) # R^2 = 0.543

df = data.frame(x,y)
p <- ggplot( df, aes( x = x/1000, y = y)) + 
  geom_point(colour = "#2b8cbe") +
  geom_smooth(colour = "#7bccc4", fill = "#ccebc5", method = 'lm') +
  scale_y_continuous(limits = c(-0.02, 0.13)) +
  #ggtitle( "rpoB, cor = 0.7365868, p-value = 3.470893e-06") +
  labs( x = "Pairwise distance (km)", y = "Fst" ) + 
  theme_bw() 

p
#ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_group_rpob.pdf', plot = p, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_group_rpobtest.pdf', plot = p, width = 7, height = 7, unit = 'cm', useDingbats = FALSE)

#=============================#
#             recA            #
#=============================#
df <- read.table("recA_input_norm.txt", header = TRUE, check.names = FALSE)   

#DKO only 
y = nrow(df)
df <- df[-4625:-y,]

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

# ======== MANTEL TEST ======# 
geodist <- as.dist(distmat)
gendist <- as.dist(fst_matrix)
mantel.rtest(geodist, gendist, nrepet = 5000) # Observation: 0.6166599, Simulated p-value: 0.01939612  
#slope 
summary(lm(gendist~geodist))


# ======== Plotting ======# 
# plotting using ggplot 
diag(distmat) = NA 
diag(fst_matrix) = NA
x <- as.vector( t( distmat ))
y <- as.vector(t(fst_matrix))
correlation(x,y) 
xx <- x/1000
mlm = lm(y ~ xx) # slope = 0.0004156
str(summary(mlm)) # R^2 = 0.38

df = data.frame(x,y)
p <- ggplot( df, aes( x = x/1000, y = y)) + 
  geom_point(colour = "#2b8cbe") +
  geom_smooth(colour = "#7bccc4", fill = "#ccebc5", method = 'lm') +
  scale_y_continuous(limits = c(-0.02, 0.13)) +
  #ggtitle( "recA, cor = 0.6166599, p-value = 0.0002845249") +
  labs( x = "Pairwise distance (km)", y = "Fst" ) + 
  theme_bw() 

p
#ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_group_reca.pdf', plot = p, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_group_recatest.pdf', plot = p, width = 7, height = 7, unit = 'cm', useDingbats = FALSE)

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

# ======== MANTEL TEST ======# 
geodist <- as.dist(distmat)
gendist <- as.dist(fst_matrix)
mantel.rtest(geodist, gendist, nrepet = 5000) # Observation: 0.6316436, Simulated p-value: 0.04659068
#slope 
summary(lm(gendist~geodist))


# ======== Plotting ======# 
# plotting using ggplot 
diag(distmat) = NA 
diag(fst_matrix) = NA
x <- as.vector( t( distmat ))
y <- as.vector(t(fst_matrix))
correlation(x,y) 
xx <- x/1000
mlm = lm(y ~ xx) # slope = 0.0001748
str(summary(mlm)) # R^2 = 0.399

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

# ======== MANTEL TEST ======# 
geodist <- as.dist(distmat)
gendist <- as.dist(fst_matrix)
mantel.rtest(geodist, gendist, nrepet = 5000) # Observation: 0.4204567, Simulated p-value: 0.114977
#slope 
summary(lm(gendist~geodist))


# ======== Plotting ======# 
# plotting using ggplot 
diag(distmat) = NA 
diag(fst_matrix) = NA
x <- as.vector( t( distmat ))
y <- as.vector(t(fst_matrix))
correlation(x,y) 
xx <- x/1000
mlm = lm(y ~ xx) # slope = 7.058e-05
str(summary(mlm)) # R^2 = 0.177

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

#DKO only 
y = nrow(df)
df <- df[-3914:-y,]

ind <- as.character(df$tree_id) # use later with adegenet (individual labels)
population <- as.character(df$sample) # use later with adegenet (population labels)

#To convert df to a genind object (adegenet), the input should only contain genotypes
x = ncol(df)
locus   <- df[, 5:x] 
df1 <- df2genind(locus, ploidy = 1, ind.names = ind, pop = population, sep = "")

# Pairwise Fst
#genet.dist(df1, method = "WC84")
fst_sample_matrix <- as.matrix(genet.dist(df1, method = "WC84"))
#fst_sample_matrix[fst_sample_matrix < 0] <- 0

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
xx <- x/1000
lm(y ~ xx) # slope = 0.0001113

df = data.frame(x,y)
p <- ggplot( df, aes( x = x/1000, y = y)) + 
  geom_point(colour = "#2b8cbe") +
  geom_smooth(colour = "#7bccc4", fill = "#ccebc5", method = 'lm') +
  scale_y_continuous(limits = c(0, 0.22)) +
  #ggtitle( "rpoB, cor = 0.1076277, p-value = 3.297469e-05") +
  labs( x = "Pairwise distance (km)", y = "Fst" ) + 
  theme_bw() 

p
#ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_rpob.pdf', plot = p, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_rpoBtest.pdf', plot = p, width = 7, height = 7, unit = 'cm', useDingbats = FALSE)


#=============================#
#             recA            #
#=============================#
df <- read.table("recA_input_norm.txt", header = TRUE, check.names = FALSE)   

#DKO only 
y = nrow(df)
df <- df[-4625:-y,]

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
xx <- x/1000
lm(y ~ xx) # slope = 1.313e-05

df = data.frame(x,y)
p <- ggplot( df, aes( x = x/1000, y = y)) + 
  geom_point(colour = "#2b8cbe") +
  geom_smooth(colour = "#7bccc4", fill = "#ccebc5", method = 'lm') +
  scale_y_continuous(limits = c(0, 0.22)) +
  #ggtitle( "recA, cor = 0.01655588, p-value = 0.4515438") +
  labs( x = "Pairwise distance (km)", y = "Fst" ) + 
  theme_bw() 
p
#ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_recA.pdf', plot = p, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/fst_dist_recAtest.pdf', plot = p, width = 7, height = 7, unit = 'cm', useDingbats = FALSE)


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
xx <- x/1000
lm(y ~ xx) # slope = 1.55e-05

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
xx <- x/1000
lm = lm(y ~ xx) # slope = -4.491e-06 
str(summary(lm))

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


