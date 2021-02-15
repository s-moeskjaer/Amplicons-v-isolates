setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/')
# packages
library(FinePop)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
#install.packages("geosphere")
library(geosphere)
library(agricolae)

#### Create data frame for analysis ####
df <- read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/ANI/ani_lastlasts.csv', sep = ';', header = FALSE)

gps = read.csv('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/ANI/DKO_GPS_ani.csv', sep = ';')
gps <- data.frame(gps$long,gps$lat) 
distmat <- distm(gps, gps, fun = distVincentyEllipsoid)


# plotting using ggplot 
diag(distmat) = NA 
diag(df) = NA
x <- as.vector(t(distmat))
y <- as.vector(t(df))
correlation(x,y) # p-value = 2.10143e-12 

df = data.frame(x,y)
p <- ggplot( df, aes( x = x, y = y)) + 
  geom_point(colour = "#4d9221") +
  ggtitle( "cor = -0.03158548, p-value = 0.004696862") +
  labs( x = "Pairwise distance (m)", y = "ANI" ) + 
  theme_bw() 

p
ggsave('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/ANI/ANI_dist_DKO.pdf', plot = p, width = 15, height = 15, unit = 'cm', useDingbats = FALSE)

