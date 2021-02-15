setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/')

#### Create data frame for analysis ####
data = read.csv('UMI_clover_order_nonodD2.csv', sep = ';')
head(data)
names(data)[names(data)=="X..."] <- "sampleID"
head(data)

#overall sample overview to add field names for DKO
envir <- read.table('zz_sample_overview_latlong.csv', header = T, sep = ";")
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

#=========================================================#
#                                                         #
#                        ALL (DKO)                        #
#                                                         #
#=========================================================#
{
#=============================#
#             rpoB            #
#=============================#
#Add annotation row based on each dataset
rpoBfilt$index <- as.numeric(row.names(rpoBfilt))
rpoBfilt <- rpoBfilt[order(rpoBfilt$index), ]
samples <- rpoBfilt$index
rpoBfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]

rpoBfilt$origin <- d2$Origin
#UK
UK <- subset(rpoBfilt, origin == 'UK')
UK_list <- colnames(UK[, colSums(UK != 0) > 0])
nonUK <- subset(rpoBfilt, origin != 'UK')
non_UK_list <- colnames(nonUK[, colSums(nonUK != 0) > 0])
x <- length(UK_list) - length(intersect(non_UK_list,UK_list))
y <- length(UK_list) -2
  
#DKO
DKO <-  subset(rpoBfilt, origin == 'DKO')
DKO_list <- colnames(DKO[, colSums(DKO != 0) > 0])
nonDKO <- subset(rpoBfilt, origin != 'DKO')
non_DKO_list <- colnames(nonDKO[, colSums(nonDKO != 0) > 0])
x <- length(DKO_list) - length(intersect(non_DKO_list,DKO_list))
y <- length(DKO_list) -2

#F
Fr <-  subset(rpoBfilt, origin == 'F')
F_list <- colnames(Fr[, colSums(Fr != 0) > 0])
nonF <- subset(rpoBfilt, origin != 'F')
non_F_list <- colnames(nonF[, colSums(nonF != 0) > 0])
x <- length(F_list) - length(intersect(non_F_list,F_list))
y <- length(F_list) -2

#DK
DK <-  subset(rpoBfilt, origin == 'DK')
DK_list <- colnames(DK[, colSums(DK != 0) > 0])
nonDK <- subset(rpoBfilt, origin != 'DK')
non_DK_list <- colnames(nonDK[, colSums(nonDK != 0) > 0])
x <- length(DK_list) - length(intersect(non_DK_list,DK_list))
y <- length(DK_list) -2

#=============================#
#             recA            #
#=============================#
#Add annotation row based on each dataset
recAfilt$index <- as.numeric(row.names(recAfilt))
recAfilt <- recAfilt[order(recAfilt$index), ]
samples <- recAfilt$index
recAfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]

recAfilt$origin <- d2$Origin
#UK
UK <- subset(recAfilt, origin == 'UK')
UK_list <- colnames(UK[, colSums(UK != 0) > 0])
nonUK <- subset(recAfilt, origin != 'UK')
non_UK_list <- colnames(nonUK[, colSums(nonUK != 0) > 0])
x <- length(UK_list) - length(intersect(non_UK_list,UK_list))
y <- length(UK_list) -2

#DKO
DKO <-  subset(recAfilt, origin == 'DKO')
DKO_list <- colnames(DKO[, colSums(DKO != 0) > 0])
nonDKO <- subset(recAfilt, origin != 'DKO')
non_DKO_list <- colnames(nonDKO[, colSums(nonDKO != 0) > 0])
x <- length(DKO_list) - length(intersect(non_DKO_list,DKO_list))
y <- length(DKO_list) -2

#F
Fr <-  subset(recAfilt, origin == 'F')
F_list <- colnames(Fr[, colSums(Fr != 0) > 0])
nonF <- subset(recAfilt, origin != 'F')
non_F_list <- colnames(nonF[, colSums(nonF != 0) > 0])
x <- length(F_list) - length(intersect(non_F_list,F_list))
y <- length(F_list) -2

#DK
DK <-  subset(recAfilt, origin == 'DK')
DK_list <- colnames(DK[, colSums(DK != 0) > 0])
nonDK <- subset(recAfilt, origin != 'DK')
non_DK_list <- colnames(nonDK[, colSums(nonDK != 0) > 0])
x <- length(DK_list) - length(intersect(non_DK_list,DK_list))
y <- length(DK_list) -2

#=============================#
#             nodA            #
#=============================#
#Add annotation row based on each dataset
nodAfilt$index <- as.numeric(row.names(nodAfilt))
nodAfilt <- nodAfilt[order(nodAfilt$index), ]
samples <- nodAfilt$index
nodAfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]

nodAfilt$origin <- d2$Origin
#UK
UK <- subset(nodAfilt, origin == 'UK')
UK_list <- colnames(UK[, colSums(UK != 0) > 0])
nonUK <- subset(nodAfilt, origin != 'UK')
non_UK_list <- colnames(nonUK[, colSums(nonUK != 0) > 0])
x <- length(UK_list) - length(intersect(non_UK_list,UK_list))
y <- length(UK_list) -2

#DKO
DKO <-  subset(nodAfilt, origin == 'DKO')
DKO_list <- colnames(DKO[, colSums(DKO != 0) > 0])
nonDKO <- subset(nodAfilt, origin != 'DKO')
non_DKO_list <- colnames(nonDKO[, colSums(nonDKO != 0) > 0])
x <- length(DKO_list) - length(intersect(non_DKO_list,DKO_list))
y <- length(DKO_list) -2

#F
Fr <-  subset(nodAfilt, origin == 'F')
F_list <- colnames(Fr[, colSums(Fr != 0) > 0])
nonF <- subset(nodAfilt, origin != 'F')
non_F_list <- colnames(nonF[, colSums(nonF != 0) > 0])
x <- length(F_list) - length(intersect(non_F_list,F_list))
y <- length(F_list) -2

#DK
DK <-  subset(nodAfilt, origin == 'DK')
DK_list <- colnames(DK[, colSums(DK != 0) > 0])
nonDK <- subset(nodAfilt, origin != 'DK')
non_DK_list <- colnames(nonDK[, colSums(nonDK != 0) > 0])
x <- length(DK_list) - length(intersect(non_DK_list,DK_list))
y <- length(DK_list) -2

#=============================#
#             nodD            #
#=============================#
#Add annotation row based on each dataset
nodDfilt$index <- as.numeric(row.names(nodDfilt))
nodDfilt <- nodDfilt[order(nodDfilt$index), ]
samples <- nodDfilt$index
nodDfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]

nodDfilt$origin <- d2$Origin
#UK
UK <- subset(nodDfilt, origin == 'UK')
UK_list <- colnames(UK[, colSums(UK != 0) > 0])
nonUK <- subset(nodDfilt, origin != 'UK')
non_UK_list <- colnames(nonUK[, colSums(nonUK != 0) > 0])
x <- length(UK_list) - length(intersect(non_UK_list,UK_list))
y <- length(UK_list) -2

#DKO
DKO <-  subset(nodDfilt, origin == 'DKO')
DKO_list <- colnames(DKO[, colSums(DKO != 0) > 0])
nonDKO <- subset(nodDfilt, origin != 'DKO')
non_DKO_list <- colnames(nonDKO[, colSums(nonDKO != 0) > 0])
x <- length(DKO_list) - length(intersect(non_DKO_list,DKO_list))
y <- length(DKO_list) -2

#F
Fr <-  subset(nodDfilt, origin == 'F')
F_list <- colnames(Fr[, colSums(Fr != 0) > 0])
nonF <- subset(nodDfilt, origin != 'F')
non_F_list <- colnames(nonF[, colSums(nonF != 0) > 0])
x <- length(F_list) - length(intersect(non_F_list,F_list))
y <- length(F_list) -2

#DK
DK <-  subset(nodDfilt, origin == 'DK')
DK_list <- colnames(DK[, colSums(DK != 0) > 0])
nonDK <- subset(nodDfilt, origin != 'DK')
non_DK_list <- colnames(nonDK[, colSums(nonDK != 0) > 0])
x <- length(DK_list) - length(intersect(non_DK_list,DK_list))
y <- length(DK_list) -2
}


#=========================================================#
#                                                         #
#                           DLF                           #
#                                                         #
#=========================================================#
{
#=============================#
#             rpoB            #
#=============================#
#Add annotation row based on each dataset
rpoBfilt$index <- as.numeric(row.names(rpoBfilt))
rpoBfilt <- rpoBfilt[order(rpoBfilt$index), ]
samples <- rpoBfilt$index
rpoBfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
rpoBfilt$origin <- d2$Origin
total_allele <- subset(rpoBfilt, origin != 'DKO')
total_allele <- colnames(total_allele[, colSums(total_allele != 0) > 0])
z <- length(total_allele) -3

#UK
UK <- subset(rpoBfilt, origin == 'UK')
UK_list <- colnames(UK[, colSums(UK != 0) > 0])
nonUK <- subset(rpoBfilt, origin != 'UK')
nonUK <- subset(nonUK, origin != 'DKO')
non_UK_list <- colnames(nonUK[, colSums(nonUK != 0) > 0])
x <- length(UK_list) - length(intersect(non_UK_list,UK_list))

#F
Fr <-  subset(rpoBfilt, origin == 'F')
F_list <- colnames(Fr[, colSums(Fr != 0) > 0])
nonF <- subset(rpoBfilt, origin != 'F')
nonF <- subset(nonF, origin != 'DKO')
non_F_list <- colnames(nonF[, colSums(nonF != 0) > 0])
x <- length(F_list) - length(intersect(non_F_list,F_list))

#DK
DK <-  subset(rpoBfilt, origin == 'DK')
DK_list <- colnames(DK[, colSums(DK != 0) > 0])
nonDK <- subset(rpoBfilt, origin != 'DK')
nonDK <- subset(nonDK, origin != 'DKO')
non_DK_list <- colnames(nonDK[, colSums(nonDK != 0) > 0])
x <- length(DK_list) - length(intersect(non_DK_list,DK_list))


#=============================#
#             recA            #
#=============================#
#Add annotation row based on each dataset
recAfilt$index <- as.numeric(row.names(recAfilt))
recAfilt <- recAfilt[order(recAfilt$index), ]
samples <- recAfilt$index
recAfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
recAfilt$origin <- d2$Origin
total_allele <- subset(recAfilt, origin != 'DKO')
total_allele <- colnames(total_allele[, colSums(total_allele != 0) > 0])
z <- length(total_allele) -3

#UK
UK <- subset(recAfilt, origin == 'UK')
UK_list <- colnames(UK[, colSums(UK != 0) > 0])
nonUK <- subset(recAfilt, origin != 'UK')
nonUK <- subset(nonUK, origin != 'DKO')
non_UK_list <- colnames(nonUK[, colSums(nonUK != 0) > 0])
x <- length(UK_list) - length(intersect(non_UK_list,UK_list))

#F
Fr <-  subset(recAfilt, origin == 'F')
F_list <- colnames(Fr[, colSums(Fr != 0) > 0])
nonF <- subset(recAfilt, origin != 'F')
nonF <- subset(nonF, origin != 'DKO')
non_F_list <- colnames(nonF[, colSums(nonF != 0) > 0])
x <- length(F_list) - length(intersect(non_F_list,F_list))

#DK
DK <-  subset(recAfilt, origin == 'DK')
DK_list <- colnames(DK[, colSums(DK != 0) > 0])
nonDK <- subset(recAfilt, origin != 'DK')
nonDK <- subset(nonDK, origin != 'DKO')
non_DK_list <- colnames(nonDK[, colSums(nonDK != 0) > 0])
x <- length(DK_list) - length(intersect(non_DK_list,DK_list))


#=============================#
#             nodA            #
#=============================#
#Add annotation row based on each dataset
nodAfilt$index <- as.numeric(row.names(nodAfilt))
nodAfilt <- nodAfilt[order(nodAfilt$index), ]
samples <- nodAfilt$index
nodAfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
nodAfilt$origin <- d2$Origin
total_allele <- subset(nodAfilt, origin != 'DKO')
total_allele <- colnames(total_allele[, colSums(total_allele != 0) > 0])
z <- length(total_allele) -3

#UK
UK <- subset(nodAfilt, origin == 'UK')
UK_list <- colnames(UK[, colSums(UK != 0) > 0])
nonUK <- subset(nodAfilt, origin != 'UK')
nonUK <- subset(nonUK, origin != 'DKO')
non_UK_list <- colnames(nonUK[, colSums(nonUK != 0) > 0])
x <- length(UK_list) - length(intersect(non_UK_list,UK_list))

#F
Fr <-  subset(nodAfilt, origin == 'F')
F_list <- colnames(Fr[, colSums(Fr != 0) > 0])
nonF <- subset(nodAfilt, origin != 'F')
nonF <- subset(nonF, origin != 'DKO')
non_F_list <- colnames(nonF[, colSums(nonF != 0) > 0])
x <- length(F_list) - length(intersect(non_F_list,F_list))

#DK
DK <-  subset(nodAfilt, origin == 'DK')
DK_list <- colnames(DK[, colSums(DK != 0) > 0])
nonDK <- subset(nodAfilt, origin != 'DK')
nonDK <- subset(nonDK, origin != 'DKO')
non_DK_list <- colnames(nonDK[, colSums(nonDK != 0) > 0])
x <- length(DK_list) - length(intersect(non_DK_list,DK_list))


#=============================#
#             nodD            #
#=============================#
#Add annotation row based on each dataset
nodDfilt$index <- as.numeric(row.names(nodDfilt))
nodDfilt <- nodDfilt[order(nodDfilt$index), ]
samples <- nodDfilt$index
nodDfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
nodDfilt$origin <- d2$Origin
total_allele <- subset(nodDfilt, origin != 'DKO')
total_allele <- colnames(total_allele[, colSums(total_allele != 0) > 0])
z <- length(total_allele) -3

#UK
UK <- subset(nodDfilt, origin == 'UK')
UK_list <- colnames(UK[, colSums(UK != 0) > 0])
nonUK <- subset(nodDfilt, origin != 'UK')
nonUK <- subset(nonUK, origin != 'DKO')
non_UK_list <- colnames(nonUK[, colSums(nonUK != 0) > 0])
x <- length(UK_list) - length(intersect(non_UK_list,UK_list))

#F
Fr <-  subset(nodDfilt, origin == 'F')
F_list <- colnames(Fr[, colSums(Fr != 0) > 0])
nonF <- subset(nodDfilt, origin != 'F')
nonF <- subset(nonF, origin != 'DKO')
non_F_list <- colnames(nonF[, colSums(nonF != 0) > 0])
x <- length(F_list) - length(intersect(non_F_list,F_list))

#DK
DK <-  subset(nodDfilt, origin == 'DK')
DK_list <- colnames(DK[, colSums(DK != 0) > 0])
nonDK <- subset(nodDfilt, origin != 'DK')
nonDK <- subset(nonDK, origin != 'DKO')
non_DK_list <- colnames(nonDK[, colSums(nonDK != 0) > 0])
x <- length(DK_list) - length(intersect(non_DK_list,DK_list))
}



#=========================================================#
#                                                         #
#              DKO groupings (incl DLF)                   #
#                                                         #
#=========================================================#
{
#=============================#
#             rpoB            #
#=============================#
#Add annotation row based on each dataset
rpoBfilt$index <- as.numeric(row.names(rpoBfilt))
rpoBfilt <- rpoBfilt[order(rpoBfilt$index), ]
samples <- rpoBfilt$index
rpoBfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
rpoBfilt$origin <- d2$Origin
rpoBfilt$grouping <- d2$Grouping

#UK
UK <- subset(rpoBfilt, origin == 'UK')
UK_list <- colnames(UK[, colSums(UK != 0) > 0])
nonUK <- subset(rpoBfilt, origin != 'UK')
non_UK_list <- colnames(nonUK[, colSums(nonUK != 0) > 0])
x <- length(UK_list) - length(intersect(non_UK_list,UK_list))

#F
Fr <-  subset(rpoBfilt, origin == 'F')
F_list <- colnames(Fr[, colSums(Fr != 0) > 0])
nonF <- subset(rpoBfilt, origin != 'F')
non_F_list <- colnames(nonF[, colSums(nonF != 0) > 0])
x <- length(F_list) - length(intersect(non_F_list,F_list))

#DK
DK <-  subset(rpoBfilt, origin == 'DK')
DK_list <- colnames(DK[, colSums(DK != 0) > 0])
nonDK <- subset(rpoBfilt, origin != 'DK')
non_DK_list <- colnames(nonDK[, colSums(nonDK != 0) > 0])
x <- length(DK_list) - length(intersect(non_DK_list,DK_list))

#DKO_1
DKO_1 <-  subset(rpoBfilt, grouping == 'DKO_1')
DKO_1_list <- colnames(DKO_1[, colSums(DKO_1 != 0) > 0])
nonDKO_1 <- subset(rpoBfilt, grouping != 'DKO_1')
non_DKO_1_list <- colnames(nonDKO_1[, colSums(nonDKO_1 != 0) > 0])
x <- length(DKO_1_list) - length(intersect(non_DKO_1_list,DKO_1_list))
y <- length(DKO_1_list) - 3

#DKO_2
DKO_2 <-  subset(rpoBfilt, grouping == 'DKO_2')
DKO_2_list <- colnames(DKO_2[, colSums(DKO_2 != 0) > 0])
nonDKO_2 <- subset(rpoBfilt, grouping != 'DKO_2')
non_DKO_2_list <- colnames(nonDKO_2[, colSums(nonDKO_2 != 0) > 0])
x <- length(DKO_2_list) - length(intersect(non_DKO_2_list,DKO_2_list))
y <- length(DKO_2_list) - 3

#DKO_3
DKO_3 <-  subset(rpoBfilt, grouping == 'DKO_3')
DKO_3_list <- colnames(DKO_3[, colSums(DKO_3 != 0) > 0])
nonDKO_3 <- subset(rpoBfilt, grouping != 'DKO_3')
non_DKO_3_list <- colnames(nonDKO_3[, colSums(nonDKO_3 != 0) > 0])
x <- length(DKO_3_list) - length(intersect(non_DKO_3_list,DKO_3_list))
y <- length(DKO_3_list) - 3

#DKO_4
DKO_4 <-  subset(rpoBfilt, grouping == 'DKO_4')
DKO_4_list <- colnames(DKO_4[, colSums(DKO_4 != 0) > 0])
nonDKO_4 <- subset(rpoBfilt, grouping != 'DKO_4')
non_DKO_4_list <- colnames(nonDKO_4[, colSums(nonDKO_4 != 0) > 0])
x <- length(DKO_4_list) - length(intersect(non_DKO_4_list,DKO_4_list))
y <- length(DKO_4_list) - 3

#DKO_5
DKO_5 <-  subset(rpoBfilt, grouping == 'DKO_5')
DKO_5_list <- colnames(DKO_5[, colSums(DKO_5 != 0) > 0])
nonDKO_5 <- subset(rpoBfilt, grouping != 'DKO_5')
non_DKO_5_list <- colnames(nonDKO_5[, colSums(nonDKO_5 != 0) > 0])
x <- length(DKO_5_list) - length(intersect(non_DKO_5_list,DKO_5_list))
y <- length(DKO_5_list) - 3 

#DKO_6
DKO_6 <-  subset(rpoBfilt, grouping == 'DKO_6')
DKO_6_list <- colnames(DKO_6[, colSums(DKO_6 != 0) > 0])
nonDKO_6 <- subset(rpoBfilt, grouping != 'DKO_6')
non_DKO_6_list <- colnames(nonDKO_6[, colSums(nonDKO_6 != 0) > 0])
x <- length(DKO_6_list) - length(intersect(non_DKO_6_list,DKO_6_list))
y <- length(DKO_6_list) - 3 



#=============================#
#             recA            #
#=============================#
#Add annotation row based on each dataset
recAfilt$index <- as.numeric(row.names(recAfilt))
recAfilt <- recAfilt[order(recAfilt$index), ]
samples <- recAfilt$index
recAfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
recAfilt$origin <- d2$Origin
recAfilt$grouping <- d2$Grouping

#DKO_1
DKO_1 <-  subset(recAfilt, grouping == 'DKO_1')
DKO_1_list <- colnames(DKO_1[, colSums(DKO_1 != 0) > 0])
nonDKO_1 <- subset(recAfilt, grouping != 'DKO_1')
non_DKO_1_list <- colnames(nonDKO_1[, colSums(nonDKO_1 != 0) > 0])
x <- length(DKO_1_list) - length(intersect(non_DKO_1_list,DKO_1_list))
y <- length(DKO_1_list) - 3

#DKO_2
DKO_2 <-  subset(recAfilt, grouping == 'DKO_2')
DKO_2_list <- colnames(DKO_2[, colSums(DKO_2 != 0) > 0])
nonDKO_2 <- subset(recAfilt, grouping != 'DKO_2')
non_DKO_2_list <- colnames(nonDKO_2[, colSums(nonDKO_2 != 0) > 0])
x <- length(DKO_2_list) - length(intersect(non_DKO_2_list,DKO_2_list))
y <- length(DKO_2_list) - 3

#DKO_3
DKO_3 <-  subset(recAfilt, grouping == 'DKO_3')
DKO_3_list <- colnames(DKO_3[, colSums(DKO_3 != 0) > 0])
nonDKO_3 <- subset(recAfilt, grouping != 'DKO_3')
non_DKO_3_list <- colnames(nonDKO_3[, colSums(nonDKO_3 != 0) > 0])
x <- length(DKO_3_list) - length(intersect(non_DKO_3_list,DKO_3_list))
y <- length(DKO_3_list) - 3

#DKO_4
DKO_4 <-  subset(recAfilt, grouping == 'DKO_4')
DKO_4_list <- colnames(DKO_4[, colSums(DKO_4 != 0) > 0])
nonDKO_4 <- subset(recAfilt, grouping != 'DKO_4')
non_DKO_4_list <- colnames(nonDKO_4[, colSums(nonDKO_4 != 0) > 0])
x <- length(DKO_4_list) - length(intersect(non_DKO_4_list,DKO_4_list))
y <- length(DKO_4_list) - 3

#DKO_5
DKO_5 <-  subset(recAfilt, grouping == 'DKO_5')
DKO_5_list <- colnames(DKO_5[, colSums(DKO_5 != 0) > 0])
nonDKO_5 <- subset(recAfilt, grouping != 'DKO_5')
non_DKO_5_list <- colnames(nonDKO_5[, colSums(nonDKO_5 != 0) > 0])
x <- length(DKO_5_list) - length(intersect(non_DKO_5_list,DKO_5_list))
y <- length(DKO_5_list) - 3 

#DKO_6
DKO_6 <-  subset(recAfilt, grouping == 'DKO_6')
DKO_6_list <- colnames(DKO_6[, colSums(DKO_6 != 0) > 0])
nonDKO_6 <- subset(recAfilt, grouping != 'DKO_6')
non_DKO_6_list <- colnames(nonDKO_6[, colSums(nonDKO_6 != 0) > 0])
x <- length(DKO_6_list) - length(intersect(non_DKO_6_list,DKO_6_list))
y <- length(DKO_6_list) - 3 




#=============================#
#             nodA            #
#=============================#
#Add annotation row based on each dataset
nodAfilt$index <- as.numeric(row.names(nodAfilt))
nodAfilt <- nodAfilt[order(nodAfilt$index), ]
samples <- nodAfilt$index
nodAfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
nodAfilt$origin <- d2$Origin
nodAfilt$grouping <- d2$Grouping

#DKO_1
DKO_1 <-  subset(nodAfilt, grouping == 'DKO_1')
DKO_1_list <- colnames(DKO_1[, colSums(DKO_1 != 0) > 0])
nonDKO_1 <- subset(nodAfilt, grouping != 'DKO_1')
non_DKO_1_list <- colnames(nonDKO_1[, colSums(nonDKO_1 != 0) > 0])
x <- length(DKO_1_list) - length(intersect(non_DKO_1_list,DKO_1_list))
y <- length(DKO_1_list) - 3

#DKO_2
DKO_2 <-  subset(nodAfilt, grouping == 'DKO_2')
DKO_2_list <- colnames(DKO_2[, colSums(DKO_2 != 0) > 0])
nonDKO_2 <- subset(nodAfilt, grouping != 'DKO_2')
non_DKO_2_list <- colnames(nonDKO_2[, colSums(nonDKO_2 != 0) > 0])
x <- length(DKO_2_list) - length(intersect(non_DKO_2_list,DKO_2_list))
y <- length(DKO_2_list) - 3

#DKO_3
DKO_3 <-  subset(nodAfilt, grouping == 'DKO_3')
DKO_3_list <- colnames(DKO_3[, colSums(DKO_3 != 0) > 0])
nonDKO_3 <- subset(nodAfilt, grouping != 'DKO_3')
non_DKO_3_list <- colnames(nonDKO_3[, colSums(nonDKO_3 != 0) > 0])
x <- length(DKO_3_list) - length(intersect(non_DKO_3_list,DKO_3_list))
y <- length(DKO_3_list) - 3

#DKO_4
DKO_4 <-  subset(nodAfilt, grouping == 'DKO_4')
DKO_4_list <- colnames(DKO_4[, colSums(DKO_4 != 0) > 0])
nonDKO_4 <- subset(nodAfilt, grouping != 'DKO_4')
non_DKO_4_list <- colnames(nonDKO_4[, colSums(nonDKO_4 != 0) > 0])
x <- length(DKO_4_list) - length(intersect(non_DKO_4_list,DKO_4_list))
y <- length(DKO_4_list) - 3

#DKO_5
DKO_5 <-  subset(nodAfilt, grouping == 'DKO_5')
DKO_5_list <- colnames(DKO_5[, colSums(DKO_5 != 0) > 0])
nonDKO_5 <- subset(nodAfilt, grouping != 'DKO_5')
non_DKO_5_list <- colnames(nonDKO_5[, colSums(nonDKO_5 != 0) > 0])
x <- length(DKO_5_list) - length(intersect(non_DKO_5_list,DKO_5_list))
y <- length(DKO_5_list) - 3 

#DKO_6
DKO_6 <-  subset(nodAfilt, grouping == 'DKO_6')
DKO_6_list <- colnames(DKO_6[, colSums(DKO_6 != 0) > 0])
nonDKO_6 <- subset(nodAfilt, grouping != 'DKO_6')
non_DKO_6_list <- colnames(nonDKO_6[, colSums(nonDKO_6 != 0) > 0])
x <- length(DKO_6_list) - length(intersect(non_DKO_6_list,DKO_6_list))
y <- length(DKO_6_list) - 3 


#=============================#
#             nodD            #
#=============================#
#Add annotation row based on each dataset
nodDfilt$index <- as.numeric(row.names(nodDfilt))
nodDfilt <- nodDfilt[order(nodDfilt$index), ]
samples <- nodDfilt$index
nodDfilt$index <- NULL

d2 <- envir[which(envir$order %in% samples),]
d2 <- d2[order(d2$order), ]
nodDfilt$origin <- d2$Origin
nodDfilt$grouping <- d2$Grouping

#DKO_1
DKO_1 <-  subset(nodDfilt, grouping == 'DKO_1')
DKO_1_list <- colnames(DKO_1[, colSums(DKO_1 != 0) > 0])
nonDKO_1 <- subset(nodDfilt, grouping != 'DKO_1')
non_DKO_1_list <- colnames(nonDKO_1[, colSums(nonDKO_1 != 0) > 0])
x <- length(DKO_1_list) - length(intersect(non_DKO_1_list,DKO_1_list))
y <- length(DKO_1_list) - 3

#DKO_2
DKO_2 <-  subset(nodDfilt, grouping == 'DKO_2')
DKO_2_list <- colnames(DKO_2[, colSums(DKO_2 != 0) > 0])
nonDKO_2 <- subset(nodDfilt, grouping != 'DKO_2')
non_DKO_2_list <- colnames(nonDKO_2[, colSums(nonDKO_2 != 0) > 0])
x <- length(DKO_2_list) - length(intersect(non_DKO_2_list,DKO_2_list))
y <- length(DKO_2_list) - 3

#DKO_3
DKO_3 <-  subset(nodDfilt, grouping == 'DKO_3')
DKO_3_list <- colnames(DKO_3[, colSums(DKO_3 != 0) > 0])
nonDKO_3 <- subset(nodDfilt, grouping != 'DKO_3')
non_DKO_3_list <- colnames(nonDKO_3[, colSums(nonDKO_3 != 0) > 0])
x <- length(DKO_3_list) - length(intersect(non_DKO_3_list,DKO_3_list))
y <- length(DKO_3_list) - 3

#DKO_4
DKO_4 <-  subset(nodDfilt, grouping == 'DKO_4')
DKO_4_list <- colnames(DKO_4[, colSums(DKO_4 != 0) > 0])
nonDKO_4 <- subset(nodDfilt, grouping != 'DKO_4')
non_DKO_4_list <- colnames(nonDKO_4[, colSums(nonDKO_4 != 0) > 0])
x <- length(DKO_4_list) - length(intersect(non_DKO_4_list,DKO_4_list))
y <- length(DKO_4_list) - 3

#DKO_5
DKO_5 <-  subset(nodDfilt, grouping == 'DKO_5')
DKO_5_list <- colnames(DKO_5[, colSums(DKO_5 != 0) > 0])
nonDKO_5 <- subset(nodDfilt, grouping != 'DKO_5')
non_DKO_5_list <- colnames(nonDKO_5[, colSums(nonDKO_5 != 0) > 0])
x <- length(DKO_5_list) - length(intersect(non_DKO_5_list,DKO_5_list))
y <- length(DKO_5_list) - 3 

#DKO_6
DKO_6 <-  subset(nodDfilt, grouping == 'DKO_6')
DKO_6_list <- colnames(DKO_6[, colSums(DKO_6 != 0) > 0])
nonDKO_6 <- subset(nodDfilt, grouping != 'DKO_6')
non_DKO_6_list <- colnames(nonDKO_6[, colSums(nonDKO_6 != 0) > 0])
x <- length(DKO_6_list) - length(intersect(non_DKO_6_list,DKO_6_list))
y <- length(DKO_6_list) - 3 
}



#=========================================================#
#                                                         #
#                  DKO groupings only                     #
#                                                         #
#=========================================================#
{
  #=============================#
  #             rpoB            #
  #=============================#
  #Add annotation row based on each dataset
  rpoBfilt$index <- as.numeric(row.names(rpoBfilt))
  rpoBfilt <- rpoBfilt[order(rpoBfilt$index), ]
  samples <- rpoBfilt$index
  rpoBfilt$index <- NULL
  
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  rpoBfilt$origin <- d2$Origin
  rpoBfilt$grouping <- d2$Grouping
  rpoBfilt1 <- subset(rpoBfilt, origin == 'DKO')
  
  #DKO_1
  DKO_1 <-  subset(rpoBfilt1, grouping == 'DKO_1')
  DKO_1_list <- colnames(DKO_1[, colSums(DKO_1 != 0) > 0])
  nonDKO_1 <- subset(rpoBfilt1, grouping != 'DKO_1')
  non_DKO_1_list <- colnames(nonDKO_1[, colSums(nonDKO_1 != 0) > 0])
  x <- length(DKO_1_list) - length(intersect(non_DKO_1_list,DKO_1_list))
  y <- length(DKO_1_list) - 3
  
  #DKO_2
  DKO_2 <-  subset(rpoBfilt1, grouping == 'DKO_2')
  DKO_2_list <- colnames(DKO_2[, colSums(DKO_2 != 0) > 0])
  nonDKO_2 <- subset(rpoBfilt1, grouping != 'DKO_2')
  non_DKO_2_list <- colnames(nonDKO_2[, colSums(nonDKO_2 != 0) > 0])
  x <- length(DKO_2_list) - length(intersect(non_DKO_2_list,DKO_2_list))
  y <- length(DKO_2_list) - 3
  
  #DKO_3
  DKO_3 <-  subset(rpoBfilt1, grouping == 'DKO_3')
  DKO_3_list <- colnames(DKO_3[, colSums(DKO_3 != 0) > 0])
  nonDKO_3 <- subset(rpoBfilt1, grouping != 'DKO_3')
  non_DKO_3_list <- colnames(nonDKO_3[, colSums(nonDKO_3 != 0) > 0])
  x <- length(DKO_3_list) - length(intersect(non_DKO_3_list,DKO_3_list))
  y <- length(DKO_3_list) - 3
  
  #DKO_4
  DKO_4 <-  subset(rpoBfilt1, grouping == 'DKO_4')
  DKO_4_list <- colnames(DKO_4[, colSums(DKO_4 != 0) > 0])
  nonDKO_4 <- subset(rpoBfilt1, grouping != 'DKO_4')
  non_DKO_4_list <- colnames(nonDKO_4[, colSums(nonDKO_4 != 0) > 0])
  x <- length(DKO_4_list) - length(intersect(non_DKO_4_list,DKO_4_list))
  y <- length(DKO_4_list) - 3
  
  #DKO_5
  DKO_5 <-  subset(rpoBfilt1, grouping == 'DKO_5')
  DKO_5_list <- colnames(DKO_5[, colSums(DKO_5 != 0) > 0])
  nonDKO_5 <- subset(rpoBfilt1, grouping != 'DKO_5')
  non_DKO_5_list <- colnames(nonDKO_5[, colSums(nonDKO_5 != 0) > 0])
  x <- length(DKO_5_list) - length(intersect(non_DKO_5_list,DKO_5_list))
  y <- length(DKO_5_list) - 3 
  
  #DKO_6
  DKO_6 <-  subset(rpoBfilt1, grouping == 'DKO_6')
  DKO_6_list <- colnames(DKO_6[, colSums(DKO_6 != 0) > 0])
  nonDKO_6 <- subset(rpoBfilt1, grouping != 'DKO_6')
  non_DKO_6_list <- colnames(nonDKO_6[, colSums(nonDKO_6 != 0) > 0])
  x <- length(DKO_6_list) - length(intersect(non_DKO_6_list,DKO_6_list))
  y <- length(DKO_6_list) - 3 
  
  
  
  #=============================#
  #             recA            #
  #=============================#
  #Add annotation row based on each dataset
  recAfilt$index <- as.numeric(row.names(recAfilt))
  recAfilt <- recAfilt[order(recAfilt$index), ]
  samples <- recAfilt$index
  recAfilt$index <- NULL
  
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  recAfilt$origin <- d2$Origin
  recAfilt$grouping <- d2$Grouping
  recAfilt1 <- subset(recAfilt, origin == 'DKO')
  
  #DKO_1
  DKO_1 <-  subset(recAfilt1, grouping == 'DKO_1')
  DKO_1_list <- colnames(DKO_1[, colSums(DKO_1 != 0) > 0])
  nonDKO_1 <- subset(recAfilt1, grouping != 'DKO_1')
  non_DKO_1_list <- colnames(nonDKO_1[, colSums(nonDKO_1 != 0) > 0])
  x <- length(DKO_1_list) - length(intersect(non_DKO_1_list,DKO_1_list))
  y <- length(DKO_1_list) - 3
  
  #DKO_2
  DKO_2 <-  subset(recAfilt1, grouping == 'DKO_2')
  DKO_2_list <- colnames(DKO_2[, colSums(DKO_2 != 0) > 0])
  nonDKO_2 <- subset(recAfilt1, grouping != 'DKO_2')
  non_DKO_2_list <- colnames(nonDKO_2[, colSums(nonDKO_2 != 0) > 0])
  x <- length(DKO_2_list) - length(intersect(non_DKO_2_list,DKO_2_list))
  y <- length(DKO_2_list) - 3
  
  #DKO_3
  DKO_3 <-  subset(recAfilt1, grouping == 'DKO_3')
  DKO_3_list <- colnames(DKO_3[, colSums(DKO_3 != 0) > 0])
  nonDKO_3 <- subset(recAfilt1, grouping != 'DKO_3')
  non_DKO_3_list <- colnames(nonDKO_3[, colSums(nonDKO_3 != 0) > 0])
  x <- length(DKO_3_list) - length(intersect(non_DKO_3_list,DKO_3_list))
  y <- length(DKO_3_list) - 3
  
  #DKO_4
  DKO_4 <-  subset(recAfilt1, grouping == 'DKO_4')
  DKO_4_list <- colnames(DKO_4[, colSums(DKO_4 != 0) > 0])
  nonDKO_4 <- subset(recAfilt1, grouping != 'DKO_4')
  non_DKO_4_list <- colnames(nonDKO_4[, colSums(nonDKO_4 != 0) > 0])
  x <- length(DKO_4_list) - length(intersect(non_DKO_4_list,DKO_4_list))
  y <- length(DKO_4_list) - 3
  
  #DKO_5
  DKO_5 <-  subset(recAfilt1, grouping == 'DKO_5')
  DKO_5_list <- colnames(DKO_5[, colSums(DKO_5 != 0) > 0])
  nonDKO_5 <- subset(recAfilt1, grouping != 'DKO_5')
  non_DKO_5_list <- colnames(nonDKO_5[, colSums(nonDKO_5 != 0) > 0])
  x <- length(DKO_5_list) - length(intersect(non_DKO_5_list,DKO_5_list))
  y <- length(DKO_5_list) - 3 
  
  #DKO_6
  DKO_6 <-  subset(recAfilt1, grouping == 'DKO_6')
  DKO_6_list <- colnames(DKO_6[, colSums(DKO_6 != 0) > 0])
  nonDKO_6 <- subset(recAfilt1, grouping != 'DKO_6')
  non_DKO_6_list <- colnames(nonDKO_6[, colSums(nonDKO_6 != 0) > 0])
  x <- length(DKO_6_list) - length(intersect(non_DKO_6_list,DKO_6_list))
  y <- length(DKO_6_list) - 3 
  
  
  
  
  #=============================#
  #             nodA            #
  #=============================#
  #Add annotation row based on each dataset
  nodAfilt$index <- as.numeric(row.names(nodAfilt))
  nodAfilt <- nodAfilt[order(nodAfilt$index), ]
  samples <- nodAfilt$index
  nodAfilt$index <- NULL
  
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  nodAfilt$origin <- d2$Origin
  nodAfilt$grouping <- d2$Grouping
  nodAfilt1 <- subset(nodAfilt, origin == 'DKO')
  
  #DKO_1
  DKO_1 <-  subset(nodAfilt1, grouping == 'DKO_1')
  DKO_1_list <- colnames(DKO_1[, colSums(DKO_1 != 0) > 0])
  nonDKO_1 <- subset(nodAfilt1, grouping != 'DKO_1')
  non_DKO_1_list <- colnames(nonDKO_1[, colSums(nonDKO_1 != 0) > 0])
  x <- length(DKO_1_list) - length(intersect(non_DKO_1_list,DKO_1_list))
  y <- length(DKO_1_list) - 3
  
  #DKO_2
  DKO_2 <-  subset(nodAfilt1, grouping == 'DKO_2')
  DKO_2_list <- colnames(DKO_2[, colSums(DKO_2 != 0) > 0])
  nonDKO_2 <- subset(nodAfilt1, grouping != 'DKO_2')
  non_DKO_2_list <- colnames(nonDKO_2[, colSums(nonDKO_2 != 0) > 0])
  x <- length(DKO_2_list) - length(intersect(non_DKO_2_list,DKO_2_list))
  y <- length(DKO_2_list) - 3
  
  #DKO_3
  DKO_3 <-  subset(nodAfilt1, grouping == 'DKO_3')
  DKO_3_list <- colnames(DKO_3[, colSums(DKO_3 != 0) > 0])
  nonDKO_3 <- subset(nodAfilt1, grouping != 'DKO_3')
  non_DKO_3_list <- colnames(nonDKO_3[, colSums(nonDKO_3 != 0) > 0])
  x <- length(DKO_3_list) - length(intersect(non_DKO_3_list,DKO_3_list))
  y <- length(DKO_3_list) - 3
  
  #DKO_4
  DKO_4 <-  subset(nodAfilt1, grouping == 'DKO_4')
  DKO_4_list <- colnames(DKO_4[, colSums(DKO_4 != 0) > 0])
  nonDKO_4 <- subset(nodAfilt1, grouping != 'DKO_4')
  non_DKO_4_list <- colnames(nonDKO_4[, colSums(nonDKO_4 != 0) > 0])
  x <- length(DKO_4_list) - length(intersect(non_DKO_4_list,DKO_4_list))
  y <- length(DKO_4_list) - 3
  
  #DKO_5
  DKO_5 <-  subset(nodAfilt1, grouping == 'DKO_5')
  DKO_5_list <- colnames(DKO_5[, colSums(DKO_5 != 0) > 0])
  nonDKO_5 <- subset(nodAfilt1, grouping != 'DKO_5')
  non_DKO_5_list <- colnames(nonDKO_5[, colSums(nonDKO_5 != 0) > 0])
  x <- length(DKO_5_list) - length(intersect(non_DKO_5_list,DKO_5_list))
  y <- length(DKO_5_list) - 3 
  
  #DKO_6
  DKO_6 <-  subset(nodAfilt1, grouping == 'DKO_6')
  DKO_6_list <- colnames(DKO_6[, colSums(DKO_6 != 0) > 0])
  nonDKO_6 <- subset(nodAfilt1, grouping != 'DKO_6')
  non_DKO_6_list <- colnames(nonDKO_6[, colSums(nonDKO_6 != 0) > 0])
  x <- length(DKO_6_list) - length(intersect(non_DKO_6_list,DKO_6_list))
  y <- length(DKO_6_list) - 3 
  
  
  #=============================#
  #             nodD            #
  #=============================#
  #Add annotation row based on each dataset
  nodDfilt$index <- as.numeric(row.names(nodDfilt))
  nodDfilt <- nodDfilt[order(nodDfilt$index), ]
  samples <- nodDfilt$index
  nodDfilt$index <- NULL
  
  d2 <- envir[which(envir$order %in% samples),]
  d2 <- d2[order(d2$order), ]
  nodDfilt$origin <- d2$Origin
  nodDfilt$grouping <- d2$Grouping
  nodDfilt1 <- subset(nodDfilt, origin == 'DKO')
  
  
  #DKO_1
  DKO_1 <-  subset(nodDfilt1, grouping == 'DKO_1')
  DKO_1_list <- colnames(DKO_1[, colSums(DKO_1 != 0) > 0])
  nonDKO_1 <- subset(nodDfilt1, grouping != 'DKO_1')
  non_DKO_1_list <- colnames(nonDKO_1[, colSums(nonDKO_1 != 0) > 0])
  x <- length(DKO_1_list) - length(intersect(non_DKO_1_list,DKO_1_list))
  y <- length(DKO_1_list) - 3
  
  #DKO_2
  DKO_2 <-  subset(nodDfilt1, grouping == 'DKO_2')
  DKO_2_list <- colnames(DKO_2[, colSums(DKO_2 != 0) > 0])
  nonDKO_2 <- subset(nodDfilt1, grouping != 'DKO_2')
  non_DKO_2_list <- colnames(nonDKO_2[, colSums(nonDKO_2 != 0) > 0])
  x <- length(DKO_2_list) - length(intersect(non_DKO_2_list,DKO_2_list))
  y <- length(DKO_2_list) - 3
  
  #DKO_3
  DKO_3 <-  subset(nodDfilt1, grouping == 'DKO_3')
  DKO_3_list <- colnames(DKO_3[, colSums(DKO_3 != 0) > 0])
  nonDKO_3 <- subset(nodDfilt1, grouping != 'DKO_3')
  non_DKO_3_list <- colnames(nonDKO_3[, colSums(nonDKO_3 != 0) > 0])
  x <- length(DKO_3_list) - length(intersect(non_DKO_3_list,DKO_3_list))
  y <- length(DKO_3_list) - 3
  
  #DKO_4
  DKO_4 <-  subset(nodDfilt1, grouping == 'DKO_4')
  DKO_4_list <- colnames(DKO_4[, colSums(DKO_4 != 0) > 0])
  nonDKO_4 <- subset(nodDfilt1, grouping != 'DKO_4')
  non_DKO_4_list <- colnames(nonDKO_4[, colSums(nonDKO_4 != 0) > 0])
  x <- length(DKO_4_list) - length(intersect(non_DKO_4_list,DKO_4_list))
  y <- length(DKO_4_list) - 3
  
  #DKO_5
  DKO_5 <-  subset(nodDfilt1, grouping == 'DKO_5')
  DKO_5_list <- colnames(DKO_5[, colSums(DKO_5 != 0) > 0])
  nonDKO_5 <- subset(nodDfilt1, grouping != 'DKO_5')
  non_DKO_5_list <- colnames(nonDKO_5[, colSums(nonDKO_5 != 0) > 0])
  x <- length(DKO_5_list) - length(intersect(non_DKO_5_list,DKO_5_list))
  y <- length(DKO_5_list) - 3 
  
  #DKO_6
  DKO_6 <-  subset(nodDfilt1, grouping == 'DKO_6')
  DKO_6_list <- colnames(DKO_6[, colSums(DKO_6 != 0) > 0])
  nonDKO_6 <- subset(nodDfilt1, grouping != 'DKO_6')
  non_DKO_6_list <- colnames(nonDKO_6[, colSums(nonDKO_6 != 0) > 0])
  x <- length(DKO_6_list) - length(intersect(non_DKO_6_list,DKO_6_list))
  y <- length(DKO_6_list) - 3 
  
}