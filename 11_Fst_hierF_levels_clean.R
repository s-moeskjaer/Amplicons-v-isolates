# https://popgen.nescent.org/DifferentiationSNP.html
# https://www.researchgate.net/publication/6074121_A_step-by-step_tutorial_to_use_HierFstat_to_analyse_populations_hierarchically_structured_at_multiple_levels
# http://www.t-de-meeus.fr/examplehier.txt

library(adegenet)
library(hierfstat)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

setwd('/Users/sara91/Documents/PhD/Materials_and_methods/QQAD/201909_MAUI_2015/Fst/fst_by_distance/groupings/hierfstat/')


#==================================================================================================================#
#                                                                                                                  #
#                                                      rpoB                                                        #
#                                                                                                                  #
#==================================================================================================================#
{
  # ========================================== #
  #            Preparing the data              #
  # ========================================== #
  {
    df <- read.table("rpoB_input_norm.txt", header = TRUE, check.names = FALSE)   
    
    # =============== LEVEL TEST =============== #
    # Add DKO/DLF level 
    df$sample[1:3913] <- "DKO"
    y = nrow(df)
    df$sample[3914:y] <- "Field trial"
    levels <- data.frame(df[, c(2:4)])
    x = ncol(df)
    loci <- df[, 5:x] 
    
    
  {loci$`1`[loci$`1`=="C"] <- 1
    loci$`1`[loci$`1`=="T"] <- 2
    loci$`2`[loci$`2`=="C"] <- 1
    loci$`2`[loci$`2`=="T"] <- 2
    loci$`3`[loci$`3`=="G"] <- 1
    loci$`3`[loci$`3`=="A"] <- 2
    loci$`4`[loci$`4`=="G"] <- 1
    loci$`4`[loci$`4`=="C"] <- 2
    loci$`5`[loci$`5`=="G"] <- 1
    loci$`5`[loci$`5`=="C"] <- 2
    loci$`6`[loci$`6`=="T"] <- 1
    loci$`6`[loci$`6`=="C"] <- 2
    loci$`7`[loci$`7`=="T"] <- 1
    loci$`7`[loci$`7`=="C"] <- 2
    loci$`8`[loci$`8`=="C"] <- 1
    loci$`8`[loci$`8`=="T"] <- 2
    loci$`9`[loci$`9`=="C"] <- 1
    loci$`9`[loci$`9`=="T"] <- 2
    loci$`10`[loci$`10`=="T"] <- 1
    loci$`10`[loci$`10`=="C"] <- 2
    loci$`11`[loci$`11`=="C"] <- 1
    loci$`11`[loci$`11`=="T"] <- 2
    loci$`12`[loci$`12`=="C"] <- 1
    loci$`12`[loci$`12`=="T"] <- 2
    loci$`12`[loci$`12`=="G"] <- 3
    loci$`13`[loci$`13`=="C"] <- 1
    loci$`13`[loci$`13`=="T"] <- 2
    loci$`14`[loci$`14`=="C"] <- 1
    loci$`14`[loci$`14`=="G"] <- 2
    loci$`15`[loci$`15`=="C"] <- 1
    loci$`15`[loci$`15`=="T"] <- 2
    loci$`16`[loci$`16`=="C"] <- 1
    loci$`16`[loci$`16`=="T"] <- 2
    loci$`17`[loci$`17`=="G"] <- 1
    loci$`17`[loci$`17`=="T"] <- 2
    loci$`18`[loci$`18`=="G"] <- 1
    loci$`18`[loci$`18`=="A"] <- 2
  }
    
    # Testing requires levels to be numerical
    levels$sample[levels$sample=="DKO"] <- 1
    levels$sample[levels$sample=="Field trial"] <- 2
    
    levels$grouping[levels$grouping=="DKO_1"] <- 1
    levels$grouping[levels$grouping=="DKO_2"] <- 2
    levels$grouping[levels$grouping=="DKO_3"] <- 3
    levels$grouping[levels$grouping=="DKO_4"] <- 4
    levels$grouping[levels$grouping=="DKO_5"] <- 5
    levels$grouping[levels$grouping=="DKO_6"] <- 6
    levels$grouping[levels$grouping=="DK"] <- 7
    levels$grouping[levels$grouping=="F"] <- 8
    levels$grouping[levels$grouping=="UK"] <- 9
    
    
    levels$grouping <- as.numeric(levels$grouping)
    levels$sample <- as.numeric(levels$sample)
    
    # ONLY FOR DLF Add clover genotype file for DLF
    {
      clover_genotype <- read.csv("clover_genotypes_DLF.csv", sep = ";")
      levels$order <- as.numeric(row.names(levels))
      levels <- merge(levels, clover_genotype, by = "sample_grouping")
      levels <- levels[order(levels$order),]
      row.names(levels) <- levels$order
      
      
      levels$Same_number_for_DLF <- as.numeric(levels$Same_number_for_DLF)
      levels$different_between_trials <- as.numeric(levels$different_between_trials)
      levels$field_number <- as.numeric(levels$field_number)
    }
    
    {loci$`1` <- as.numeric(loci$`1`)
    loci$`2` <- as.numeric(loci$`2`)
    loci$`3` <- as.numeric(loci$`3`)
    loci$`4` <- as.numeric(loci$`4`)
    loci$`5` <- as.numeric(loci$`5`)
    loci$`6` <- as.numeric(loci$`6`)
    loci$`7` <- as.numeric(loci$`7`)
    loci$`8` <- as.numeric(loci$`8`)
    loci$`9` <- as.numeric(loci$`9`)
    loci$`10` <- as.numeric(loci$`10`)
    loci$`11` <- as.numeric(loci$`11`)
    loci$`12` <- as.numeric(loci$`12`)
    loci$`13` <- as.numeric(loci$`13`)
    loci$`14` <- as.numeric(loci$`14`)
    loci$`15` <- as.numeric(loci$`15`)
    loci$`16` <- as.numeric(loci$`16`)
    loci$`17` <- as.numeric(loci$`17`)
    loci$`18` <- as.numeric(loci$`18`)
    }
    
    {levels$sample_grouping[levels$sample_grouping=="DKO_1_DKO_1"] <- 1
    levels$sample_grouping[levels$sample_grouping=="DKO_2_DKO_1"] <- 2
    levels$sample_grouping[levels$sample_grouping=="DKO_4_DKO_1"] <- 3
    levels$sample_grouping[levels$sample_grouping=="DKO_5_DKO_1"] <- 4
    levels$sample_grouping[levels$sample_grouping=="DKO_6_DKO_1"] <- 5
    levels$sample_grouping[levels$sample_grouping=="DKO_7_DKO_1"] <- 6
    levels$sample_grouping[levels$sample_grouping=="DKO_8_DKO_1"] <- 7
    levels$sample_grouping[levels$sample_grouping=="DKO_9_DKO_1"] <- 8
    levels$sample_grouping[levels$sample_grouping=="DKO_10_DKO_1"] <- 9
    levels$sample_grouping[levels$sample_grouping=="DKO_11_DKO_1"] <- 10
    levels$sample_grouping[levels$sample_grouping=="DKO_12_DKO_1"] <- 11
    levels$sample_grouping[levels$sample_grouping=="DKO_13_DKO_1"] <- 12
    levels$sample_grouping[levels$sample_grouping=="DKO_14_DKO_1"] <- 13
    levels$sample_grouping[levels$sample_grouping=="DKO_15_DKO_4"] <- 14  
    levels$sample_grouping[levels$sample_grouping=="DKO_16_DKO_4"] <- 15
    levels$sample_grouping[levels$sample_grouping=="DKO_17_DKO_4"] <- 16
    levels$sample_grouping[levels$sample_grouping=="DKO_18_DKO_4"] <- 17
    levels$sample_grouping[levels$sample_grouping=="DKO_19_DKO_4"] <- 18
    levels$sample_grouping[levels$sample_grouping=="DKO_21_DKO_4"] <- 19
    levels$sample_grouping[levels$sample_grouping=="DKO_22_DKO_4"] <- 20
    levels$sample_grouping[levels$sample_grouping=="DKO_23_DKO_4"] <- 21
    levels$sample_grouping[levels$sample_grouping=="DKO_25_DKO_5"] <- 22
    levels$sample_grouping[levels$sample_grouping=="DKO_26_DKO_5"] <- 23
    levels$sample_grouping[levels$sample_grouping=="DKO_27_DKO_3"] <- 24
    levels$sample_grouping[levels$sample_grouping=="DKO_28_DKO_4"] <- 25
    levels$sample_grouping[levels$sample_grouping=="DKO_29_DKO_4"] <- 26
    levels$sample_grouping[levels$sample_grouping=="DKO_31_DKO_4"] <- 27
    levels$sample_grouping[levels$sample_grouping=="DKO_35_DKO_6"] <- 28
    levels$sample_grouping[levels$sample_grouping=="DKO_36_DKO_6"] <- 29
    levels$sample_grouping[levels$sample_grouping=="DKO_37_DKO_2"] <- 30
    levels$sample_grouping[levels$sample_grouping=="DKO_38_DKO_2"] <- 31
    levels$sample_grouping[levels$sample_grouping=="DKO_39_DKO_2"] <- 32
    levels$sample_grouping[levels$sample_grouping=="DKO_40_DKO_5"] <- 33
    levels$sample_grouping[levels$sample_grouping=="DKO_41_DKO_3"] <- 34
    levels$sample_grouping[levels$sample_grouping=="DKO_42_DKO_3"] <- 35
    levels$sample_grouping[levels$sample_grouping=="DKO_43_DKO_5"] <- 36
    levels$sample_grouping[levels$sample_grouping=="DKO_44_DKO_4"] <- 37
    levels$sample_grouping[levels$sample_grouping=="DKO_49_DKO_2"] <- 38
    levels$sample_grouping[levels$sample_grouping=="DKO_50_DKO_6"] <- 39
    levels$sample_grouping[levels$sample_grouping=="DK_131_DK"] <- 40
    levels$sample_grouping[levels$sample_grouping=="DK_132_DK"] <- 41
    levels$sample_grouping[levels$sample_grouping=="DK_133_DK"] <- 42
    levels$sample_grouping[levels$sample_grouping=="DK_134_DK"] <- 43
    levels$sample_grouping[levels$sample_grouping=="DK_135_DK"] <- 44
    levels$sample_grouping[levels$sample_grouping=="DK_136_DK"] <- 45
    levels$sample_grouping[levels$sample_grouping=="DK_137_DK"] <- 46
    levels$sample_grouping[levels$sample_grouping=="DK_138_DK"] <- 47
    levels$sample_grouping[levels$sample_grouping=="DK_162_DK"] <- 48
    levels$sample_grouping[levels$sample_grouping=="DK_163_DK"] <- 49
    levels$sample_grouping[levels$sample_grouping=="DK_164_DK"] <- 50
    levels$sample_grouping[levels$sample_grouping=="DK_165_DK"] <- 51
    levels$sample_grouping[levels$sample_grouping=="DK_166_DK"] <- 52
    levels$sample_grouping[levels$sample_grouping=="DK_169_DK"] <- 53 
    levels$sample_grouping[levels$sample_grouping=="DK_170_DK"] <- 54
    levels$sample_grouping[levels$sample_grouping=="F_51_F"] <- 55
    levels$sample_grouping[levels$sample_grouping=="F_52_F"] <- 56
    levels$sample_grouping[levels$sample_grouping=="F_53_F"] <- 57
    levels$sample_grouping[levels$sample_grouping=="F_54_F"] <- 58
    levels$sample_grouping[levels$sample_grouping=="F_55_F"] <- 59
    levels$sample_grouping[levels$sample_grouping=="F_56_F"] <- 60
    levels$sample_grouping[levels$sample_grouping=="F_57_F"] <- 61
    levels$sample_grouping[levels$sample_grouping=="F_58_F"] <- 62
    levels$sample_grouping[levels$sample_grouping=="F_59_F"] <- 63
    levels$sample_grouping[levels$sample_grouping=="F_60_F"] <- 64
    levels$sample_grouping[levels$sample_grouping=="F_61_F"] <- 65
    levels$sample_grouping[levels$sample_grouping=="F_62_F"] <- 66
    levels$sample_grouping[levels$sample_grouping=="F_63_F"] <- 67
    levels$sample_grouping[levels$sample_grouping=="F_64_F"] <- 68
    levels$sample_grouping[levels$sample_grouping=="F_65_F"] <- 69
    levels$sample_grouping[levels$sample_grouping=="F_66_F"] <- 70
    levels$sample_grouping[levels$sample_grouping=="F_67_F"] <- 71
    levels$sample_grouping[levels$sample_grouping=="F_68_F"] <- 72
    levels$sample_grouping[levels$sample_grouping=="F_69_F"] <- 73
    levels$sample_grouping[levels$sample_grouping=="F_70_F"] <- 74
    levels$sample_grouping[levels$sample_grouping=="F_71_F"] <- 75
    levels$sample_grouping[levels$sample_grouping=="F_72_F"] <- 76
    levels$sample_grouping[levels$sample_grouping=="F_73_F"] <- 78
    levels$sample_grouping[levels$sample_grouping=="F_74_F"] <- 79
    levels$sample_grouping[levels$sample_grouping=="F_75_F"] <- 80
    levels$sample_grouping[levels$sample_grouping=="F_79_F"] <- 81
    levels$sample_grouping[levels$sample_grouping=="F_80_F"] <- 82
    levels$sample_grouping[levels$sample_grouping=="F_81_F"] <- 83
    levels$sample_grouping[levels$sample_grouping=="F_82_F"] <- 84
    levels$sample_grouping[levels$sample_grouping=="F_84_F"] <- 85
    levels$sample_grouping[levels$sample_grouping=="F_85_F"] <- 86
    levels$sample_grouping[levels$sample_grouping=="F_86_F"] <- 87
    levels$sample_grouping[levels$sample_grouping=="F_87_F"] <- 88
    levels$sample_grouping[levels$sample_grouping=="F_88_F"] <- 89
    levels$sample_grouping[levels$sample_grouping=="F_90_F"] <- 90
    levels$sample_grouping[levels$sample_grouping=="UK_116_UK"] <- 91
    levels$sample_grouping[levels$sample_grouping=="UK_117_UK"] <- 92
    levels$sample_grouping[levels$sample_grouping=="UK_118_UK"] <- 93
    levels$sample_grouping[levels$sample_grouping=="UK_119_UK"] <- 94
    levels$sample_grouping[levels$sample_grouping=="UK_120_UK"] <- 95
    levels$sample_grouping[levels$sample_grouping=="UK_122_UK"] <- 96
    levels$sample_grouping[levels$sample_grouping=="UK_125_UK"] <- 97
    levels$sample_grouping[levels$sample_grouping=="UK_126_UK"] <- 98
    levels$sample_grouping[levels$sample_grouping=="UK_127_UK"] <- 99
    levels$sample_grouping[levels$sample_grouping=="UK_128_UK"] <- 100
    levels$sample_grouping[levels$sample_grouping=="UK_129_UK"] <- 101
    levels$sample_grouping[levels$sample_grouping=="UK_130_UK"] <- 102
    levels$sample_grouping[levels$sample_grouping=="UK_91_UK"] <- 103
    levels$sample_grouping[levels$sample_grouping=="UK_92_UK"] <- 104
    levels$sample_grouping[levels$sample_grouping=="UK_93_UK"] <- 105
    levels$sample_grouping[levels$sample_grouping=="UK_95_UK"] <- 106
    }
    levels$sample_grouping <- as.numeric(levels$sample_grouping)
    
  }
  
  # ========================================== #
  #                All samples                 #
  # ========================================== #
  {
    #Calculate Fst values 
    varcomp.glob(levels,loci, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of management/total (DKOvDLF)
    test.between(loci,rand.unit=levels$grouping,test=levels$sample,nperm=1000) # p-value = 0.106
    # Test for the effect of grouping/total 
    test.between(loci,rand.unit=levels$sample_grouping,test=levels$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    levels$individual_sample <- row.names(levels)
    test.between(loci,rand.unit=levels$individual_sample,test=levels$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of groupings/management, issue the command
    test.within(loci, test=levels$grouping, within = levels$sample, nperm=1000) # p-value = 0.001
    # To test for the effect of field_plot/management, issue the command
    test.within(loci, test=levels$sample_grouping, within = levels$sample, nperm=1000) # p-value = 
    # To test for the effect of field_plot/grouping, issue the command
    test.within(loci, test=levels$sample_grouping, within = levels$grouping, nperm=1000) # p-value = 0.001
    
  }
  
  # ========================================== #
  #               Field trials                 #
  # ========================================== #
  {
    #Choose field trial samples
    y = nrow(loci)
    loci2 <- loci[3914:y,]
    y = nrow(levels)
    levels2 <- levels[3914:y,]
    
    levels2 <- data.frame(levels2[, c(2:3)])
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$sample_grouping, within = levels2$grouping, nperm=1000) # p-value = 0.001
  }
  # ========================================== #
  #                  Organic                   #
  # ========================================== #
  {
    #Choose DKO samples
    y = nrow(loci)
    loci2 <- loci[-3914:-y,]
    y = nrow(levels)
    levels2 <- levels[-3914:-y,]
    
    levels2 <- data.frame(levels2[, c(2:3)])
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total 
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.004
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$sample_grouping, within = levels2$grouping, nperm=1000) # p-value = 0.001
    
  }  
  
  # ========================================== #
  #      Field trials with clover genotype     #
  # ========================================== #
  {
    #Choose field trial samples
    y = nrow(df)
    loci2 <- loci[3914:y,]
    levels2 <- levels[3914:y,]
    
    #Same number for clover genotypes between field trial sites
    levels2 <- data.frame(levels2$grouping, levels2$field_number, levels2$Same_number_for_DLF, levels2$sample_grouping)
    
    #Different number for clover genotypes between field trial sites
    #levels2 <- data.frame(levels2$grouping, levels2$field_number, levels2$different_between_trials, levels2$sample_grouping)
    
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$levels2.sample_grouping, within = levels2$levels2.grouping, nperm=1000) # p-value = 0.001
  }
}    


#==================================================================================================================#
#                                                                                                                  #
#                                                      recA                                                        #
#                                                                                                                  #
#==================================================================================================================#
{
  # ========================================== #
  #            Preparing the data              #
  # ========================================== #
  {
    df <- read.table("recA_input_norm.txt", header = TRUE, check.names = FALSE)   
    
    # =============== LEVEL TEST =============== #
    # Add DKO/DLF level 
    df$sample[1:4624] <- "DKO"
    y = nrow(df)
    df$sample[4625:y] <- "Field trial"
    levels <- data.frame(df[, c(2:4)])
    x = ncol(df)
    loci <- df[, 5:x]  
    
    {loci$`1`[loci$`1`=="C"] <- 1
    loci$`1`[loci$`1`=="T"] <- 2
    loci$`2`[loci$`2`=="C"] <- 1
    loci$`2`[loci$`2`=="T"] <- 2
    loci$`3`[loci$`3`=="C"] <- 1
    loci$`3`[loci$`3`=="T"] <- 2
    loci$`4`[loci$`4`=="C"] <- 1
    loci$`4`[loci$`4`=="T"] <- 2
    loci$`5`[loci$`5`=="A"] <- 1
    loci$`5`[loci$`5`=="G"] <- 2
    loci$`6`[loci$`6`=="A"] <- 1
    loci$`6`[loci$`6`=="G"] <- 2
    loci$`7`[loci$`7`=="C"] <- 1
    loci$`7`[loci$`7`=="T"] <- 2
    loci$`7`[loci$`7`=="A"] <- 3
    loci$`8`[loci$`8`=="C"] <- 1
    loci$`8`[loci$`8`=="T"] <- 2
    loci$`9`[loci$`9`=="A"] <- 1
    loci$`9`[loci$`9`=="G"] <- 2
    loci$`10`[loci$`10`=="T"] <- 1
    loci$`10`[loci$`10`=="C"] <- 2
    loci$`11`[loci$`11`=="C"] <- 1
    loci$`11`[loci$`11`=="T"] <- 2
    loci$`12`[loci$`12`=="A"] <- 1
    loci$`12`[loci$`12`=="G"] <- 2
    loci$`12`[loci$`12`=="G"] <- 3
    loci$`13`[loci$`13`=="C"] <- 1
    loci$`13`[loci$`13`=="T"] <- 2
    loci$`14`[loci$`14`=="T"] <- 1
    loci$`14`[loci$`14`=="C"] <- 2
    loci$`15`[loci$`15`=="C"] <- 1
    loci$`15`[loci$`15`=="T"] <- 2
    loci$`16`[loci$`16`=="C"] <- 1
    loci$`16`[loci$`16`=="T"] <- 2
    loci$`17` <- NULL
  }
    # ONLY FOR DLF Add clover genotype file for DLF
    {
      clover_genotype <- read.csv("clover_genotypes_DLF.csv", sep = ";")
      levels$order <- as.numeric(row.names(levels))
      levels <- merge(levels, clover_genotype, by = "sample_grouping")
      levels <- levels[order(levels$order),]
      row.names(levels) <- levels$order
      
      
      levels$Same_number_for_DLF <- as.numeric(levels$Same_number_for_DLF)
      levels$different_between_trials <- as.numeric(levels$different_between_trials)
      levels$field_number <- as.numeric(levels$field_number)
    }
    
    # Testing requires levels to be numerical
    {levels$sample[levels$sample=="DKO"] <- 1
    levels$sample[levels$sample=="Field trial"] <- 2
    
    levels$grouping[levels$grouping=="DKO_1"] <- 1
    levels$grouping[levels$grouping=="DKO_2"] <- 2
    levels$grouping[levels$grouping=="DKO_3"] <- 3
    levels$grouping[levels$grouping=="DKO_4"] <- 4
    levels$grouping[levels$grouping=="DKO_5"] <- 5
    levels$grouping[levels$grouping=="DKO_6"] <- 6
    levels$grouping[levels$grouping=="DK"] <- 7
    levels$grouping[levels$grouping=="F"] <- 8
    levels$grouping[levels$grouping=="UK"] <- 9
    
    
    levels$grouping <- as.numeric(levels$grouping)
    levels$sample <- as.numeric(levels$sample)
    }
    
    {loci$`1` <- as.numeric(loci$`1`)
    loci$`2` <- as.numeric(loci$`2`)
    loci$`3` <- as.numeric(loci$`3`)
    loci$`4` <- as.numeric(loci$`4`)
    loci$`5` <- as.numeric(loci$`5`)
    loci$`6` <- as.numeric(loci$`6`)
    loci$`7` <- as.numeric(loci$`7`)
    loci$`8` <- as.numeric(loci$`8`)
    loci$`9` <- as.numeric(loci$`9`)
    loci$`10` <- as.numeric(loci$`10`)
    loci$`11` <- as.numeric(loci$`11`)
    loci$`12` <- as.numeric(loci$`12`)
    loci$`13` <- as.numeric(loci$`13`)
    loci$`14` <- as.numeric(loci$`14`)
    loci$`15` <- as.numeric(loci$`15`)
    loci$`16` <- as.numeric(loci$`16`)
    }
    
    {levels$sample_grouping[levels$sample_grouping=="DKO_1_DKO_1"] <-1
    levels$sample_grouping[levels$sample_grouping=="DKO_2_DKO_1"] <-2
    levels$sample_grouping[levels$sample_grouping=="DKO_3_DKO_1"] <-3
    levels$sample_grouping[levels$sample_grouping=="DKO_4_DKO_1"] <-4
    levels$sample_grouping[levels$sample_grouping=="DKO_5_DKO_1"] <-5
    levels$sample_grouping[levels$sample_grouping=="DKO_6_DKO_1"] <-6
    levels$sample_grouping[levels$sample_grouping=="DKO_7_DKO_1"] <-7
    levels$sample_grouping[levels$sample_grouping=="DKO_8_DKO_1"] <-8
    levels$sample_grouping[levels$sample_grouping=="DKO_9_DKO_1"] <-9
    levels$sample_grouping[levels$sample_grouping=="DKO_10_DKO_1"] <-10
    levels$sample_grouping[levels$sample_grouping=="DKO_11_DKO_1"] <-11
    levels$sample_grouping[levels$sample_grouping=="DKO_12_DKO_1"] <-12
    levels$sample_grouping[levels$sample_grouping=="DKO_13_DKO_1"] <-13
    levels$sample_grouping[levels$sample_grouping=="DKO_14_DKO_1"] <-14
    levels$sample_grouping[levels$sample_grouping=="DKO_15_DKO_4"] <-15
    levels$sample_grouping[levels$sample_grouping=="DKO_16_DKO_4"] <-16
    levels$sample_grouping[levels$sample_grouping=="DKO_17_DKO_4"] <-17
    levels$sample_grouping[levels$sample_grouping=="DKO_18_DKO_4"] <-18
    levels$sample_grouping[levels$sample_grouping=="DKO_19_DKO_4"] <-19
    levels$sample_grouping[levels$sample_grouping=="DKO_20_DKO_4"] <-20
    levels$sample_grouping[levels$sample_grouping=="DKO_21_DKO_4"] <-21
    levels$sample_grouping[levels$sample_grouping=="DKO_22_DKO_4"] <-22
    levels$sample_grouping[levels$sample_grouping=="DKO_23_DKO_4"] <-23
    levels$sample_grouping[levels$sample_grouping=="DKO_24_DKO_5"] <-24
    levels$sample_grouping[levels$sample_grouping=="DKO_25_DKO_5"] <-25
    levels$sample_grouping[levels$sample_grouping=="DKO_26_DKO_5"] <-26
    levels$sample_grouping[levels$sample_grouping=="DKO_27_DKO_3"] <-27
    levels$sample_grouping[levels$sample_grouping=="DKO_28_DKO_4"] <-28
    levels$sample_grouping[levels$sample_grouping=="DKO_29_DKO_4"] <-29
    levels$sample_grouping[levels$sample_grouping=="DKO_30_DKO_4"] <-30
    levels$sample_grouping[levels$sample_grouping=="DKO_31_DKO_4"] <-31
    levels$sample_grouping[levels$sample_grouping=="DKO_32_DKO_4"] <-32
    levels$sample_grouping[levels$sample_grouping=="DKO_33_DKO_6"] <-33
    levels$sample_grouping[levels$sample_grouping=="DKO_34_DKO_6"] <-34
    levels$sample_grouping[levels$sample_grouping=="DKO_35_DKO_6"] <-35
    levels$sample_grouping[levels$sample_grouping=="DKO_36_DKO_6"] <-36
    levels$sample_grouping[levels$sample_grouping=="DKO_37_DKO_2"] <-37
    levels$sample_grouping[levels$sample_grouping=="DKO_38_DKO_2"] <-38
    levels$sample_grouping[levels$sample_grouping=="DKO_39_DKO_2"] <-39
    levels$sample_grouping[levels$sample_grouping=="DKO_40_DKO_5"] <-40
    levels$sample_grouping[levels$sample_grouping=="DKO_41_DKO_3"] <-41
    levels$sample_grouping[levels$sample_grouping=="DKO_42_DKO_3"] <-42
    levels$sample_grouping[levels$sample_grouping=="DKO_43_DKO_5"] <-43
    levels$sample_grouping[levels$sample_grouping=="DKO_44_DKO_4"] <-44
    levels$sample_grouping[levels$sample_grouping=="DKO_49_DKO_2"] <-45
    levels$sample_grouping[levels$sample_grouping=="DKO_50_DKO_6"] <-46
    levels$sample_grouping[levels$sample_grouping=="DK_131_DK"] <-47
    levels$sample_grouping[levels$sample_grouping=="DK_132_DK"] <-48
    levels$sample_grouping[levels$sample_grouping=="DK_133_DK"] <-49
    levels$sample_grouping[levels$sample_grouping=="DK_134_DK"] <-50
    levels$sample_grouping[levels$sample_grouping=="DK_135_DK"] <-51
    levels$sample_grouping[levels$sample_grouping=="DK_136_DK"] <-52
    levels$sample_grouping[levels$sample_grouping=="DK_137_DK"] <-53
    levels$sample_grouping[levels$sample_grouping=="DK_139_DK"] <-54
    levels$sample_grouping[levels$sample_grouping=="DK_140_DK"] <-55
    levels$sample_grouping[levels$sample_grouping=="DK_141_DK"] <-56
    levels$sample_grouping[levels$sample_grouping=="DK_142_DK"] <-57
    levels$sample_grouping[levels$sample_grouping=="DK_143_DK"] <-58
    levels$sample_grouping[levels$sample_grouping=="DK_144_DK"] <-59
    levels$sample_grouping[levels$sample_grouping=="DK_145_DK"] <-60
    levels$sample_grouping[levels$sample_grouping=="DK_146_DK"] <-61
    levels$sample_grouping[levels$sample_grouping=="DK_147_DK"] <-62
    levels$sample_grouping[levels$sample_grouping=="DK_148_DK"] <-63
    levels$sample_grouping[levels$sample_grouping=="DK_149_DK"] <-64
    levels$sample_grouping[levels$sample_grouping=="DK_150_DK"] <-65
    levels$sample_grouping[levels$sample_grouping=="DK_151_DK"] <-66
    levels$sample_grouping[levels$sample_grouping=="DK_152_DK"] <-67
    levels$sample_grouping[levels$sample_grouping=="DK_153_DK"] <-68
    levels$sample_grouping[levels$sample_grouping=="DK_154_DK"] <-69
    levels$sample_grouping[levels$sample_grouping=="DK_155_DK"] <-70
    levels$sample_grouping[levels$sample_grouping=="DK_156_DK"] <-71
    levels$sample_grouping[levels$sample_grouping=="DK_157_DK"] <-72
    levels$sample_grouping[levels$sample_grouping=="DK_158_DK"] <-73
    levels$sample_grouping[levels$sample_grouping=="DK_159_DK"] <-74
    levels$sample_grouping[levels$sample_grouping=="DK_160_DK"] <-75
    levels$sample_grouping[levels$sample_grouping=="DK_161_DK"] <-76
    levels$sample_grouping[levels$sample_grouping=="DK_162_DK"] <-77
    levels$sample_grouping[levels$sample_grouping=="DK_163_DK"] <-78
    levels$sample_grouping[levels$sample_grouping=="DK_164_DK"] <-79
    levels$sample_grouping[levels$sample_grouping=="DK_165_DK"] <-80
    levels$sample_grouping[levels$sample_grouping=="DK_166_DK"] <-81
    levels$sample_grouping[levels$sample_grouping=="DK_167_DK"] <-82
    levels$sample_grouping[levels$sample_grouping=="DK_169_DK"] <-83
    levels$sample_grouping[levels$sample_grouping=="DK_170_DK"] <-84
    levels$sample_grouping[levels$sample_grouping=="F_51_F"] <-85
    levels$sample_grouping[levels$sample_grouping=="F_52_F"] <-86
    levels$sample_grouping[levels$sample_grouping=="F_53_F"] <-87
    levels$sample_grouping[levels$sample_grouping=="F_54_F"] <-88
    levels$sample_grouping[levels$sample_grouping=="F_55_F"] <-89
    levels$sample_grouping[levels$sample_grouping=="F_56_F"] <-90
    levels$sample_grouping[levels$sample_grouping=="F_57_F"] <-91
    levels$sample_grouping[levels$sample_grouping=="F_58_F"] <-92
    levels$sample_grouping[levels$sample_grouping=="F_59_F"] <-93
    levels$sample_grouping[levels$sample_grouping=="F_60_F"] <-94
    levels$sample_grouping[levels$sample_grouping=="F_61_F"] <-95
    levels$sample_grouping[levels$sample_grouping=="F_62_F"] <-96
    levels$sample_grouping[levels$sample_grouping=="F_63_F"] <-97
    levels$sample_grouping[levels$sample_grouping=="F_64_F"] <-98
    levels$sample_grouping[levels$sample_grouping=="F_65_F"] <-99
    levels$sample_grouping[levels$sample_grouping=="F_66_F"] <-100
    levels$sample_grouping[levels$sample_grouping=="F_67_F"] <-101
    levels$sample_grouping[levels$sample_grouping=="F_68_F"] <-102
    levels$sample_grouping[levels$sample_grouping=="F_69_F"] <-103
    levels$sample_grouping[levels$sample_grouping=="F_70_F"] <-104
    levels$sample_grouping[levels$sample_grouping=="F_71_F"] <-105
    levels$sample_grouping[levels$sample_grouping=="F_72_F"] <-106
    levels$sample_grouping[levels$sample_grouping=="F_73_F"] <-107
    levels$sample_grouping[levels$sample_grouping=="F_74_F"] <-108
    levels$sample_grouping[levels$sample_grouping=="F_79_F"] <-109
    levels$sample_grouping[levels$sample_grouping=="F_80_F"] <-110
    levels$sample_grouping[levels$sample_grouping=="F_81_F"] <-111
    levels$sample_grouping[levels$sample_grouping=="F_82_F"] <-112
    levels$sample_grouping[levels$sample_grouping=="F_85_F"] <-113
    levels$sample_grouping[levels$sample_grouping=="F_86_F"] <-114
    levels$sample_grouping[levels$sample_grouping=="F_87_F"] <-115
    levels$sample_grouping[levels$sample_grouping=="F_88_F"] <-116
    levels$sample_grouping[levels$sample_grouping=="F_90_F"] <-117
    levels$sample_grouping[levels$sample_grouping=="UK_101_UK"] <-118
    levels$sample_grouping[levels$sample_grouping=="UK_102_UK"] <-119
    levels$sample_grouping[levels$sample_grouping=="UK_103_UK"] <-120
    levels$sample_grouping[levels$sample_grouping=="UK_104_UK"] <-121
    levels$sample_grouping[levels$sample_grouping=="UK_105_UK"] <-122
    levels$sample_grouping[levels$sample_grouping=="UK_106_UK"] <-123
    levels$sample_grouping[levels$sample_grouping=="UK_107_UK"] <-124
    levels$sample_grouping[levels$sample_grouping=="UK_108_UK"] <-125
    levels$sample_grouping[levels$sample_grouping=="UK_110_UK"] <-126
    levels$sample_grouping[levels$sample_grouping=="UK_111_UK"] <-127
    levels$sample_grouping[levels$sample_grouping=="UK_112_UK"] <-128
    levels$sample_grouping[levels$sample_grouping=="UK_113_UK"] <-129
    levels$sample_grouping[levels$sample_grouping=="UK_114_UK"] <-130
    levels$sample_grouping[levels$sample_grouping=="UK_115_UK"] <-131
    levels$sample_grouping[levels$sample_grouping=="UK_116_UK"] <-132
    levels$sample_grouping[levels$sample_grouping=="UK_117_UK"] <-133
    levels$sample_grouping[levels$sample_grouping=="UK_118_UK"] <-134
    levels$sample_grouping[levels$sample_grouping=="UK_119_UK"] <-135
    levels$sample_grouping[levels$sample_grouping=="UK_120_UK"] <-136
    levels$sample_grouping[levels$sample_grouping=="UK_121_UK"] <-137
    levels$sample_grouping[levels$sample_grouping=="UK_122_UK"] <-138
    levels$sample_grouping[levels$sample_grouping=="UK_125_UK"] <-139
    levels$sample_grouping[levels$sample_grouping=="UK_126_UK"] <-140
    levels$sample_grouping[levels$sample_grouping=="UK_127_UK"] <-141
    levels$sample_grouping[levels$sample_grouping=="UK_128_UK"] <-142
    levels$sample_grouping[levels$sample_grouping=="UK_129_UK"] <-143
    levels$sample_grouping[levels$sample_grouping=="UK_130_UK"] <-144
    levels$sample_grouping[levels$sample_grouping=="UK_91_UK"] <-145
    levels$sample_grouping[levels$sample_grouping=="UK_92_UK"] <-146
    levels$sample_grouping[levels$sample_grouping=="UK_93_UK"] <-147
    levels$sample_grouping[levels$sample_grouping=="UK_94_UK"] <-148
    levels$sample_grouping[levels$sample_grouping=="UK_95_UK"] <-149
    levels$sample_grouping[levels$sample_grouping=="UK_96_UK"] <-150
    levels$sample_grouping[levels$sample_grouping=="UK_97_UK"] <-151
    levels$sample_grouping[levels$sample_grouping=="UK_98_UK"] <-152
    levels$sample_grouping[levels$sample_grouping=="UK_99_UK"] <-153
    }
    levels$sample_grouping <- as.numeric(levels$sample_grouping)
    
  }
  
  # ========================================== #
  #                All samples                 #
  # ========================================== #
  {
    #Calculate Fst values 
    varcomp.glob(levels,loci, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of management/total (DKOvDLF)
    test.between(loci,rand.unit=levels$grouping,test=levels$sample,nperm=1000) # p-value = 0.207
    # Test for the effect of grouping/total 
    test.between(loci,rand.unit=levels$sample_grouping,test=levels$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    levels$individual_sample <- row.names(levels)
    test.between(loci,rand.unit=levels$individual_sample,test=levels$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of groupings/management, issue the command
    test.within(loci, test=levels$grouping, within = levels$sample, nperm=1000) # p-value = 0.001
    # To test for the effect of field_plot/management, issue the command
    test.within(loci, test=levels$sample_grouping, within = levels$sample, nperm=1000) # p-value = 0.001
    # To test for the effect of field_plot/grouping, issue the command
    test.within(loci, test=levels$sample_grouping, within = levels$grouping, nperm=1000) # p-value = 0.001
    
  }
  
  # ========================================== #
  #               Field trials                 #
  # ========================================== #
  {
    #Choose field trial samples
    y = nrow(df)
    loci2 <- loci[4625:y,]
    y = nrow(levels)
    levels2 <- levels[4625:y,]
    
    levels2 <- data.frame(levels2[, c(2:3)])
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$sample_grouping, within = levels2$grouping, nperm=1000) # p-value = 0.001
  }
  # ========================================== #
  #                  Organic                   #
  # ========================================== #
  {
    #Choose DKO samples
    y = nrow(df)
    loci2 <- loci[-4625:-y,]
    y = nrow(levels)
    levels2 <- levels[-4625:-y,]
    
    levels2 <- data.frame(levels2[, c(2:3)])
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total 
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.002
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.005
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$sample_grouping, within = levels2$grouping, nperm=1000) # p-value = 0.001
    
  }  
  
  # ========================================== #
  #      Field trials with clover genotype     #
  # ========================================== #
  {
    #Choose field trial samples
    y = nrow(df)
    loci2 <- loci[4625:y,]
    levels2 <- levels[4625:y,]
    
    #Same number for clover genotypes between field trial sites
    levels2 <- data.frame(levels2$grouping, levels2$field_number, levels2$Same_number_for_DLF, levels2$sample_grouping)
    
    #Different number for clover genotypes between field trial sites
    #levels2 <- data.frame(levels2$grouping, levels2$field_number, levels2$different_between_trials, levels2$sample_grouping)
    
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$levels2.sample_grouping, within = levels2$levels2.grouping, nperm=1000) # p-value = 0.001
  }
}  



#==================================================================================================================#
#                                                                                                                  #
#                                                      nodA                                                        #
#                                                                                                                  #
#==================================================================================================================#
{
  # ========================================== #
  #            Preparing the data              #
  # ========================================== #
  {
    df <- read.table("nodA_input_norm.txt", header = TRUE, check.names = FALSE)   
    
    # =============== LEVEL TEST =============== #
    # Add DKO/DLF level 
    df$sample[1:3413] <- "DKO"
    y = nrow(df)
    df$sample[3414:y] <- "Field trial"
    levels <- data.frame(df[, c(2:4)])
    x = ncol(df)
    
    loci <- df[, 5:x] 
  
    {loci$`1`[loci$`1`=="C"] <- 1
    loci$`1`[loci$`1`=="A"] <- 2
    loci$`2` <-NULL
    loci$`3` <-NULL
    loci$`4`[loci$`4`=="T"] <- 1
    loci$`4`[loci$`4`=="G"] <- 2
    loci$`5`[loci$`5`=="A"] <- 1
    loci$`5`[loci$`5`=="G"] <- 2
    loci$`6` <-NULL
    loci$`7`[loci$`7`=="C"] <- 1
    loci$`7`[loci$`7`=="T"] <- 2
    loci$`8`[loci$`8`=="T"] <- 1
    loci$`8`[loci$`8`=="C"] <- 2
    loci$`9`[loci$`9`=="T"] <- 1
    loci$`9`[loci$`9`=="C"] <- 2
    loci$`10`[loci$`10`=="G"] <- 1
    loci$`10`[loci$`10`=="A"] <- 2
    loci$`11`[loci$`11`=="T"] <- 1
    loci$`11`[loci$`11`=="C"] <- 2
    loci$`12`[loci$`12`=="T"] <- 1
    loci$`12`[loci$`12`=="C"] <- 2
    loci$`13`[loci$`13`=="G"] <- 1
    loci$`13`[loci$`13`=="T"] <- 2
    loci$`14`[loci$`14`=="G"] <- 1
    loci$`14`[loci$`14`=="A"] <- 2
    loci$`15`[loci$`15`=="C"] <- 1
    loci$`15`[loci$`15`=="A"] <- 2
    loci$`16`[loci$`16`=="C"] <- 1
    loci$`16`[loci$`16`=="G"] <- 2
    loci$`17`[loci$`17`=="T"] <- 1
    loci$`17`[loci$`17`=="C"] <- 2
    loci$`18`[loci$`18`=="C"] <- 1
    loci$`18`[loci$`18`=="G"] <- 2
    loci$`19`[loci$`19`=="G"] <- 1
    loci$`19`[loci$`19`=="A"] <- 2
    loci$`20`[loci$`20`=="C"] <- 1
    loci$`20`[loci$`20`=="T"] <- 2
    loci$`21`[loci$`21`=="G"] <- 1
    loci$`21`[loci$`21`=="A"] <- 2
    loci$`22`[loci$`22`=="G"] <- 1
    loci$`22`[loci$`22`=="A"] <- 2
    loci$`23`[loci$`23`=="G"] <- 1
    loci$`23`[loci$`23`=="A"] <- 2
    loci$`24`[loci$`24`=="T"] <- 1
    loci$`24`[loci$`24`=="C"] <- 2
    loci$`25`[loci$`25`=="A"] <- 1
    loci$`25`[loci$`25`=="G"] <- 2
  }
    
    # ONLY FOR DLF Add clover genotype file for DLF
    {
      clover_genotype <- read.csv("clover_genotypes_DLF.csv", sep = ";")
      levels$order <- as.numeric(row.names(levels))
      levels <- merge(levels, clover_genotype, by = "sample_grouping")
      levels <- levels[order(levels$order),]
      row.names(levels) <- levels$order
      
      
      levels$Same_number_for_DLF <- as.numeric(levels$Same_number_for_DLF)
      levels$different_between_trials <- as.numeric(levels$different_between_trials)
      levels$field_number <- as.numeric(levels$field_number)
    }
    # Testing requires levels to be numerical
    {levels$sample[levels$sample=="DKO"] <- 1
    levels$sample[levels$sample=="Field trial"] <- 2
    
    levels$grouping[levels$grouping=="DKO_1"] <- 1
    levels$grouping[levels$grouping=="DKO_2"] <- 2
    levels$grouping[levels$grouping=="DKO_3"] <- 3
    levels$grouping[levels$grouping=="DKO_4"] <- 4
    levels$grouping[levels$grouping=="DKO_5"] <- 5
    levels$grouping[levels$grouping=="DKO_6"] <- 6
    levels$grouping[levels$grouping=="DK"] <- 7
    levels$grouping[levels$grouping=="F"] <- 8
    levels$grouping[levels$grouping=="UK"] <- 9
    
    
    levels$grouping <- as.numeric(levels$grouping)
    levels$sample <- as.numeric(levels$sample)
    }
    
    {loci$`1` <- as.numeric(loci$`1`)
    loci$`4` <- as.numeric(loci$`4`)
    loci$`5` <- as.numeric(loci$`5`)
    loci$`7` <- as.numeric(loci$`7`)
    loci$`8` <- as.numeric(loci$`8`)
    loci$`9` <- as.numeric(loci$`9`)
    loci$`10` <- as.numeric(loci$`10`)
    loci$`11` <- as.numeric(loci$`11`)
    loci$`12` <- as.numeric(loci$`12`)
    loci$`13` <- as.numeric(loci$`13`)
    loci$`14` <- as.numeric(loci$`14`)
    loci$`15` <- as.numeric(loci$`15`)
    loci$`16` <- as.numeric(loci$`16`)
    loci$`17` <- as.numeric(loci$`17`)
    loci$`18` <- as.numeric(loci$`18`)
    loci$`19` <- as.numeric(loci$`19`)
    loci$`20` <- as.numeric(loci$`20`)
    loci$`21` <- as.numeric(loci$`21`)
    loci$`22` <- as.numeric(loci$`22`)
    loci$`23` <- as.numeric(loci$`23`)
    loci$`24` <- as.numeric(loci$`24`)
    loci$`25` <- as.numeric(loci$`25`)
    }
    
    {levels$sample_grouping[levels$sample_grouping=="DKO_1_DKO_1"] <-1
    levels$sample_grouping[levels$sample_grouping=="DKO_2_DKO_1"] <-2
    levels$sample_grouping[levels$sample_grouping=="DKO_3_DKO_1"] <-3
    levels$sample_grouping[levels$sample_grouping=="DKO_4_DKO_1"] <-4
    levels$sample_grouping[levels$sample_grouping=="DKO_5_DKO_1"] <-5
    levels$sample_grouping[levels$sample_grouping=="DKO_6_DKO_1"] <-6
    levels$sample_grouping[levels$sample_grouping=="DKO_7_DKO_1"] <-7
    levels$sample_grouping[levels$sample_grouping=="DKO_8_DKO_1"] <-8
    levels$sample_grouping[levels$sample_grouping=="DKO_9_DKO_1"] <-9
    levels$sample_grouping[levels$sample_grouping=="DKO_10_DKO_1"] <-10
    levels$sample_grouping[levels$sample_grouping=="DKO_11_DKO_1"] <-11
    levels$sample_grouping[levels$sample_grouping=="DKO_12_DKO_1"] <-12
    levels$sample_grouping[levels$sample_grouping=="DKO_13_DKO_1"] <-13
    levels$sample_grouping[levels$sample_grouping=="DKO_14_DKO_1"] <-14
    levels$sample_grouping[levels$sample_grouping=="DKO_15_DKO_4"] <-15
    levels$sample_grouping[levels$sample_grouping=="DKO_16_DKO_4"] <-16
    levels$sample_grouping[levels$sample_grouping=="DKO_17_DKO_4"] <-17
    levels$sample_grouping[levels$sample_grouping=="DKO_18_DKO_4"] <-18
    levels$sample_grouping[levels$sample_grouping=="DKO_19_DKO_4"] <-19
    levels$sample_grouping[levels$sample_grouping=="DKO_20_DKO_4"] <-20
    levels$sample_grouping[levels$sample_grouping=="DKO_21_DKO_4"] <-21
    levels$sample_grouping[levels$sample_grouping=="DKO_22_DKO_4"] <-22
    levels$sample_grouping[levels$sample_grouping=="DKO_23_DKO_4"] <-23
    levels$sample_grouping[levels$sample_grouping=="DKO_24_DKO_5"] <-24
    levels$sample_grouping[levels$sample_grouping=="DKO_25_DKO_5"] <-25
    levels$sample_grouping[levels$sample_grouping=="DKO_26_DKO_5"] <-26
    levels$sample_grouping[levels$sample_grouping=="DKO_27_DKO_3"] <-27
    levels$sample_grouping[levels$sample_grouping=="DKO_28_DKO_4"] <-28
    levels$sample_grouping[levels$sample_grouping=="DKO_29_DKO_4"] <-29
    levels$sample_grouping[levels$sample_grouping=="DKO_30_DKO_4"] <-30
    levels$sample_grouping[levels$sample_grouping=="DKO_31_DKO_4"] <-31
    levels$sample_grouping[levels$sample_grouping=="DKO_32_DKO_4"] <-32
    levels$sample_grouping[levels$sample_grouping=="DKO_33_DKO_6"] <-33
    levels$sample_grouping[levels$sample_grouping=="DKO_34_DKO_6"] <-34
    levels$sample_grouping[levels$sample_grouping=="DKO_35_DKO_6"] <-35
    levels$sample_grouping[levels$sample_grouping=="DKO_36_DKO_6"] <-36
    levels$sample_grouping[levels$sample_grouping=="DKO_37_DKO_2"] <-37
    levels$sample_grouping[levels$sample_grouping=="DKO_38_DKO_2"] <-38
    levels$sample_grouping[levels$sample_grouping=="DKO_39_DKO_2"] <-39
    levels$sample_grouping[levels$sample_grouping=="DKO_40_DKO_5"] <-40
    levels$sample_grouping[levels$sample_grouping=="DKO_41_DKO_3"] <-41
    levels$sample_grouping[levels$sample_grouping=="DKO_42_DKO_3"] <-42
    levels$sample_grouping[levels$sample_grouping=="DKO_43_DKO_5"] <-43
    levels$sample_grouping[levels$sample_grouping=="DKO_44_DKO_4"] <-44
    levels$sample_grouping[levels$sample_grouping=="DKO_49_DKO_2"] <-45
    levels$sample_grouping[levels$sample_grouping=="DKO_50_DKO_6"] <-46
    levels$sample_grouping[levels$sample_grouping=="DK_131_DK"] <-47
    levels$sample_grouping[levels$sample_grouping=="DK_132_DK"] <-48
    levels$sample_grouping[levels$sample_grouping=="DK_133_DK"] <-49
    levels$sample_grouping[levels$sample_grouping=="DK_134_DK"] <-50
    levels$sample_grouping[levels$sample_grouping=="DK_135_DK"] <-51
    levels$sample_grouping[levels$sample_grouping=="DK_136_DK"] <-52
    levels$sample_grouping[levels$sample_grouping=="DK_137_DK"] <-53
    levels$sample_grouping[levels$sample_grouping=="DK_139_DK"] <-54
    levels$sample_grouping[levels$sample_grouping=="DK_140_DK"] <-55
    levels$sample_grouping[levels$sample_grouping=="DK_141_DK"] <-56
    levels$sample_grouping[levels$sample_grouping=="DK_142_DK"] <-57
    levels$sample_grouping[levels$sample_grouping=="DK_143_DK"] <-58
    levels$sample_grouping[levels$sample_grouping=="DK_144_DK"] <-59
    levels$sample_grouping[levels$sample_grouping=="DK_145_DK"] <-60
    levels$sample_grouping[levels$sample_grouping=="DK_146_DK"] <-61
    levels$sample_grouping[levels$sample_grouping=="DK_147_DK"] <-62
    levels$sample_grouping[levels$sample_grouping=="DK_148_DK"] <-63
    levels$sample_grouping[levels$sample_grouping=="DK_149_DK"] <-64
    levels$sample_grouping[levels$sample_grouping=="DK_150_DK"] <-65
    levels$sample_grouping[levels$sample_grouping=="DK_151_DK"] <-66
    levels$sample_grouping[levels$sample_grouping=="DK_152_DK"] <-67
    levels$sample_grouping[levels$sample_grouping=="DK_153_DK"] <-68
    levels$sample_grouping[levels$sample_grouping=="DK_154_DK"] <-69
    levels$sample_grouping[levels$sample_grouping=="DK_155_DK"] <-70
    levels$sample_grouping[levels$sample_grouping=="DK_156_DK"] <-71
    levels$sample_grouping[levels$sample_grouping=="DK_157_DK"] <-72
    levels$sample_grouping[levels$sample_grouping=="DK_158_DK"] <-73
    levels$sample_grouping[levels$sample_grouping=="DK_159_DK"] <-74
    levels$sample_grouping[levels$sample_grouping=="DK_160_DK"] <-75
    levels$sample_grouping[levels$sample_grouping=="DK_161_DK"] <-76
    levels$sample_grouping[levels$sample_grouping=="DK_162_DK"] <-77
    levels$sample_grouping[levels$sample_grouping=="DK_163_DK"] <-78
    levels$sample_grouping[levels$sample_grouping=="DK_164_DK"] <-79
    levels$sample_grouping[levels$sample_grouping=="DK_165_DK"] <-80
    levels$sample_grouping[levels$sample_grouping=="DK_166_DK"] <-81
    levels$sample_grouping[levels$sample_grouping=="DK_167_DK"] <-82
    levels$sample_grouping[levels$sample_grouping=="DK_169_DK"] <-83
    levels$sample_grouping[levels$sample_grouping=="DK_170_DK"] <-84
    levels$sample_grouping[levels$sample_grouping=="F_51_F"] <-85
    levels$sample_grouping[levels$sample_grouping=="F_52_F"] <-86
    levels$sample_grouping[levels$sample_grouping=="F_53_F"] <-87
    levels$sample_grouping[levels$sample_grouping=="F_54_F"] <-88
    levels$sample_grouping[levels$sample_grouping=="F_55_F"] <-89
    levels$sample_grouping[levels$sample_grouping=="F_56_F"] <-90
    levels$sample_grouping[levels$sample_grouping=="F_57_F"] <-91
    levels$sample_grouping[levels$sample_grouping=="F_58_F"] <-92
    levels$sample_grouping[levels$sample_grouping=="F_59_F"] <-93
    levels$sample_grouping[levels$sample_grouping=="F_60_F"] <-94
    levels$sample_grouping[levels$sample_grouping=="F_61_F"] <-95
    levels$sample_grouping[levels$sample_grouping=="F_62_F"] <-96
    levels$sample_grouping[levels$sample_grouping=="F_63_F"] <-97
    levels$sample_grouping[levels$sample_grouping=="F_64_F"] <-98
    levels$sample_grouping[levels$sample_grouping=="F_65_F"] <-99
    levels$sample_grouping[levels$sample_grouping=="F_66_F"] <-100
    levels$sample_grouping[levels$sample_grouping=="F_67_F"] <-101
    levels$sample_grouping[levels$sample_grouping=="F_68_F"] <-102
    levels$sample_grouping[levels$sample_grouping=="F_69_F"] <-103
    levels$sample_grouping[levels$sample_grouping=="F_70_F"] <-104
    levels$sample_grouping[levels$sample_grouping=="F_71_F"] <-105
    levels$sample_grouping[levels$sample_grouping=="F_72_F"] <-106
    levels$sample_grouping[levels$sample_grouping=="F_73_F"] <-107
    levels$sample_grouping[levels$sample_grouping=="F_74_F"] <-108
    levels$sample_grouping[levels$sample_grouping=="F_79_F"] <-109
    levels$sample_grouping[levels$sample_grouping=="F_80_F"] <-110
    levels$sample_grouping[levels$sample_grouping=="F_81_F"] <-111
    levels$sample_grouping[levels$sample_grouping=="F_82_F"] <-112
    levels$sample_grouping[levels$sample_grouping=="F_85_F"] <-113
    levels$sample_grouping[levels$sample_grouping=="F_86_F"] <-114
    levels$sample_grouping[levels$sample_grouping=="F_87_F"] <-115
    levels$sample_grouping[levels$sample_grouping=="F_88_F"] <-116
    levels$sample_grouping[levels$sample_grouping=="F_90_F"] <-117
    levels$sample_grouping[levels$sample_grouping=="UK_101_UK"] <-118
    levels$sample_grouping[levels$sample_grouping=="UK_102_UK"] <-119
    levels$sample_grouping[levels$sample_grouping=="UK_103_UK"] <-120
    levels$sample_grouping[levels$sample_grouping=="UK_104_UK"] <-121
    levels$sample_grouping[levels$sample_grouping=="UK_105_UK"] <-122
    levels$sample_grouping[levels$sample_grouping=="UK_106_UK"] <-123
    levels$sample_grouping[levels$sample_grouping=="UK_107_UK"] <-124
    levels$sample_grouping[levels$sample_grouping=="UK_108_UK"] <-125
    levels$sample_grouping[levels$sample_grouping=="UK_109_UK"] <-126
    levels$sample_grouping[levels$sample_grouping=="UK_110_UK"] <-127
    levels$sample_grouping[levels$sample_grouping=="UK_111_UK"] <-128
    levels$sample_grouping[levels$sample_grouping=="UK_112_UK"] <-129
    levels$sample_grouping[levels$sample_grouping=="UK_113_UK"] <-130
    levels$sample_grouping[levels$sample_grouping=="UK_114_UK"] <-131
    levels$sample_grouping[levels$sample_grouping=="UK_115_UK"] <-132
    levels$sample_grouping[levels$sample_grouping=="UK_116_UK"] <-133
    levels$sample_grouping[levels$sample_grouping=="UK_117_UK"] <-134
    levels$sample_grouping[levels$sample_grouping=="UK_118_UK"] <-135
    levels$sample_grouping[levels$sample_grouping=="UK_119_UK"] <-136
    levels$sample_grouping[levels$sample_grouping=="UK_120_UK"] <-137
    levels$sample_grouping[levels$sample_grouping=="UK_121_UK"] <-138
    levels$sample_grouping[levels$sample_grouping=="UK_122_UK"] <-139
    levels$sample_grouping[levels$sample_grouping=="UK_125_UK"] <-140
    levels$sample_grouping[levels$sample_grouping=="UK_126_UK"] <-141
    levels$sample_grouping[levels$sample_grouping=="UK_127_UK"] <-142
    levels$sample_grouping[levels$sample_grouping=="UK_128_UK"] <-143
    levels$sample_grouping[levels$sample_grouping=="UK_129_UK"] <-144
    levels$sample_grouping[levels$sample_grouping=="UK_130_UK"] <-145
    levels$sample_grouping[levels$sample_grouping=="UK_91_UK"] <-146
    levels$sample_grouping[levels$sample_grouping=="UK_92_UK"] <-147
    levels$sample_grouping[levels$sample_grouping=="UK_93_UK"] <-148
    levels$sample_grouping[levels$sample_grouping=="UK_94_UK"] <-149
    levels$sample_grouping[levels$sample_grouping=="UK_95_UK"] <-150
    levels$sample_grouping[levels$sample_grouping=="UK_96_UK"] <-151
    levels$sample_grouping[levels$sample_grouping=="UK_97_UK"] <-152
    levels$sample_grouping[levels$sample_grouping=="UK_98_UK"] <-153
    levels$sample_grouping[levels$sample_grouping=="UK_99_UK"] <-154
    }
    levels$sample_grouping <- as.numeric(levels$sample_grouping)
    
  }
  
  # ========================================== #
  #                All samples                 #
  # ========================================== #
  {
    #Calculate Fst values 
    varcomp.glob(levels,loci, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of management/total (DKOvDLF)
    test.between(loci,rand.unit=levels$grouping,test=levels$sample,nperm=1000) # p-value = 0.188
    # Test for the effect of grouping/total 
    test.between(loci,rand.unit=levels$sample_grouping,test=levels$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    levels$individual_sample <- row.names(levels)
    test.between(loci,rand.unit=levels$individual_sample,test=levels$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of groupings/management, issue the command
    test.within(loci, test=levels$grouping, within = levels$sample, nperm=1000) # p-value = 0.001
    # To test for the effect of field_plot/management, issue the command
    test.within(loci, test=levels$sample_grouping, within = levels$sample, nperm=1000) # p-value = 0.001
    # To test for the effect of field_plot/grouping, issue the command
    test.within(loci, test=levels$sample_grouping, within = levels$grouping, nperm=1000) # p-value = 0.001
    
  }
  
  # ========================================== #
  #               Field trials                 #
  # ========================================== #
  {
    #Choose field trial samples
    y = nrow(df)
    loci2 <- loci[3414:y,]
    y = nrow(levels)
    levels2 <- levels[3414:y,]
    
    levels2 <- data.frame(levels2[, c(2:3)])
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$sample_grouping, within = levels2$grouping, nperm=1000) # p-value = 0.001
  }
  # ========================================== #
  #                  Organic                   #
  # ========================================== #
  {
    #Choose DKO samples
    y = nrow(df)
    loci2 <- loci[-3414:-y,]
    y = nrow(levels)
    levels2 <- levels[-3414:-y,]
    
    levels2 <- data.frame(levels2[, c(2:3)])
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total 
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.221
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.033
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$sample_grouping, within = levels2$grouping, nperm=1000) # p-value = 0.001
    
  }  
  
  # ========================================== #
  #      Field trials with clover genotype     #
  # ========================================== #
  {
    #Choose field trial samples
    y = nrow(df)
    loci2 <- loci[3414:y,]
    levels2 <- levels[3414:y,]
    
    #Same number for clover genotypes between field trial sites
    levels2 <- data.frame(levels2$grouping, levels2$field_number, levels2$Same_number_for_DLF, levels2$sample_grouping)
    
    #Different number for clover genotypes between field trial sites
    #levels2 <- data.frame(levels2$grouping, levels2$field_number, levels2$different_between_trials, levels2$sample_grouping)
    
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$levels2.sample_grouping, within = levels2$levels2.grouping, nperm=1000) # p-value = 0.001
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$levels2.field_number, within = levels2$levels2.grouping, nperm=1000) # p-value = 0.001
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$levels2.Same_number_for_DLF, within = levels2$levels2.grouping, nperm=1000) # p-value = 0.001
  }
} 


#==================================================================================================================#
#                                                                                                                  #
#                                                      nodD                                                        #
#                                                                                                                  #
#==================================================================================================================#
{
  # ========================================== #
  #            Preparing the data              #
  # ========================================== #
  {
    df <- read.table("nodD_input_norm.txt", header = TRUE, check.names = FALSE)   
    
    # =============== LEVEL TEST =============== #
    # Add DKO/DLF level 
    df$sample[1:3699] <- "DKO"
    y = nrow(df)
    df$sample[3700:y] <- "Field trial"
    levels <- data.frame(df[, c(2:4)])
    
    # Make loci numerical
    {
      x = ncol(df)
    loci <- df[, 5:x] 
    loci$`1`[loci$`1`=="T"] <- 1
    loci$`1`[loci$`1`=="C"] <- 2
    loci$`2`[loci$`2`=="G"] <- 1
    loci$`2`[loci$`2`=="C"] <- 2
    loci$`3`[loci$`3`=="G"] <- 1
    loci$`3`[loci$`3`=="T"] <- 2
    loci$`4`[loci$`4`=="G"] <- 1
    loci$`4`[loci$`4`=="A"] <- 2
    loci$`5` <-NULL
    loci$`6` <-NULL
    loci$`7`[loci$`7`=="C"] <- 1
    loci$`7`[loci$`7`=="A"] <- 2
    loci$`8`[loci$`8`=="G"] <- 1
    loci$`8`[loci$`8`=="A"] <- 2
    loci$`9`[loci$`9`=="G"] <- 1
    loci$`9`[loci$`9`=="C"] <- 2
    loci$`10`[loci$`10`=="A"] <- 1
    loci$`10`[loci$`10`=="G"] <- 2
    loci$`11`[loci$`11`=="T"] <- 1
    loci$`11`[loci$`11`=="C"] <- 2
    loci$`12`[loci$`12`=="C"] <- 1
    loci$`12`[loci$`12`=="T"] <- 2
    loci$`13`[loci$`13`=="C"] <- 1
    loci$`13`[loci$`13`=="T"] <- 2
    loci$`14`[loci$`14`=="A"] <- 1
    loci$`14`[loci$`14`=="G"] <- 2
    loci$`15`[loci$`15`=="A"] <- 1
    loci$`15`[loci$`15`=="G"] <- 2
    loci$`16`[loci$`16`=="T"] <- 1
    loci$`16`[loci$`16`=="C"] <- 2
    loci$`17`[loci$`17`=="G"] <- 1
    loci$`17`[loci$`17`=="C"] <- 2
    loci$`18`[loci$`18`=="C"] <- 1
    loci$`18`[loci$`18`=="T"] <- 2
    loci$`19`[loci$`19`=="G"] <- 1
    loci$`19`[loci$`19`=="A"] <- 2
    loci$`20`[loci$`20`=="T"] <- 1
    loci$`20`[loci$`20`=="C"] <- 2
    loci$`21`[loci$`21`=="G"] <- 1
    loci$`21`[loci$`21`=="T"] <- 2
    loci$`22` <- NULL
    loci$`23`[loci$`23`=="C"] <- 1
    loci$`23`[loci$`23`=="A"] <- 2
    loci$`24`[loci$`24`=="C"] <- 1
    loci$`24`[loci$`24`=="T"] <- 2
    loci$`25`[loci$`25`=="G"] <- 1
    loci$`25`[loci$`25`=="A"] <- 2
    loci$`26`[loci$`26`=="A"] <- 1
    loci$`26`[loci$`26`=="C"] <- 2
    loci$`27`[loci$`27`=="G"] <- 1
    loci$`27`[loci$`27`=="A"] <- 2
    loci$`28`[loci$`28`=="A"] <- 1
    loci$`28`[loci$`28`=="G"] <- 2
    loci$`29`[loci$`29`=="T"] <- 1
    loci$`29`[loci$`29`=="C"] <- 2
    loci$`30`[loci$`30`=="A"] <- 1
    loci$`30`[loci$`30`=="G"] <- 2
    loci$`31`[loci$`31`=="T"] <- 1
    loci$`31`[loci$`31`=="G"] <- 2
    loci$`32` <- NULL
    loci$`33`[loci$`33`=="T"] <- 1
    loci$`33`[loci$`33`=="A"] <- 2
    loci$`34`[loci$`34`=="G"] <- 1
    loci$`34`[loci$`34`=="T"] <- 2
    loci$`35`[loci$`35`=="G"] <- 1
    loci$`35`[loci$`35`=="C"] <- 2
    loci$`36` <- NULL
    loci$`37` <- NULL
    loci$`38`[loci$`38`=="G"] <- 1
    loci$`38`[loci$`38`=="A"] <- 2
    
    loci$`1` <- as.numeric(loci$`1`)
    loci$`2` <- as.numeric(loci$`2`)
    loci$`3` <- as.numeric(loci$`3`)
    loci$`4` <- as.numeric(loci$`4`)
    loci$`7` <- as.numeric(loci$`7`)
    loci$`8` <- as.numeric(loci$`8`)
    loci$`9` <- as.numeric(loci$`9`)
    loci$`10` <- as.numeric(loci$`10`)
    loci$`11` <- as.numeric(loci$`11`)
    loci$`12` <- as.numeric(loci$`12`)
    loci$`13` <- as.numeric(loci$`13`)
    loci$`14` <- as.numeric(loci$`14`)
    loci$`15` <- as.numeric(loci$`15`)
    loci$`16` <- as.numeric(loci$`16`)
    loci$`17` <- as.numeric(loci$`17`)
    loci$`18` <- as.numeric(loci$`18`)
    loci$`19` <- as.numeric(loci$`19`)
    loci$`20` <- as.numeric(loci$`20`)
    loci$`21` <- as.numeric(loci$`21`)
    loci$`23` <- as.numeric(loci$`23`)
    loci$`24` <- as.numeric(loci$`24`)
    loci$`25` <- as.numeric(loci$`25`)
    loci$`26` <- as.numeric(loci$`26`)
    loci$`27` <- as.numeric(loci$`27`)
    loci$`28` <- as.numeric(loci$`28`)
    loci$`29` <- as.numeric(loci$`29`)
    loci$`30` <- as.numeric(loci$`30`)
    loci$`31` <- as.numeric(loci$`31`)
    loci$`33` <- as.numeric(loci$`33`)
    loci$`34` <- as.numeric(loci$`34`)
    loci$`35` <- as.numeric(loci$`35`)
    loci$`38` <- as.numeric(loci$`38`)
    }
    
    # Testing requires levels to be numerical
    levels$sample[levels$sample=="DKO"] <- 1
    levels$sample[levels$sample=="Field trial"] <- 2
    
    levels$grouping[levels$grouping=="DKO_1"] <- 1
    levels$grouping[levels$grouping=="DKO_2"] <- 2
    levels$grouping[levels$grouping=="DKO_3"] <- 3
    levels$grouping[levels$grouping=="DKO_4"] <- 4
    levels$grouping[levels$grouping=="DKO_5"] <- 5
    levels$grouping[levels$grouping=="DKO_6"] <- 6
    levels$grouping[levels$grouping=="DK"] <- 7
    levels$grouping[levels$grouping=="F"] <- 8
    levels$grouping[levels$grouping=="UK"] <- 9
    
    
    levels$grouping <- as.numeric(levels$grouping)
    levels$sample <- as.numeric(levels$sample)
    
    # ONLY FOR DLF Add clover genotype file for DLF
    {
      clover_genotype <- read.csv("clover_genotypes_DLF.csv", sep = ";")
      levels$order <- as.numeric(row.names(levels))
      levels <- merge(levels, clover_genotype, by = "sample_grouping")
      levels <- levels[order(levels$order),]
      row.names(levels) <- levels$order
      
      
      levels$Same_number_for_DLF <- as.numeric(levels$Same_number_for_DLF)
      levels$different_between_trials <- as.numeric(levels$different_between_trials)
      levels$field_number <- as.numeric(levels$field_number)
    }
    
    # Make levels numerical
    {
    levels$sample_grouping[levels$sample_grouping=="DKO_1_DKO_1"] <- 1
    levels$sample_grouping[levels$sample_grouping=="DKO_2_DKO_1"] <- 2
    levels$sample_grouping[levels$sample_grouping=="DKO_3_DKO_1"] <- 3
    levels$sample_grouping[levels$sample_grouping=="DKO_4_DKO_1"] <- 4
    levels$sample_grouping[levels$sample_grouping=="DKO_5_DKO_1"] <- 5
    levels$sample_grouping[levels$sample_grouping=="DKO_6_DKO_1"] <- 6
    levels$sample_grouping[levels$sample_grouping=="DKO_7_DKO_1"] <- 7
    levels$sample_grouping[levels$sample_grouping=="DKO_8_DKO_1"] <- 8
    levels$sample_grouping[levels$sample_grouping=="DKO_9_DKO_1"] <- 9
    levels$sample_grouping[levels$sample_grouping=="DKO_10_DKO_1"] <- 10
    levels$sample_grouping[levels$sample_grouping=="DKO_11_DKO_1"] <- 11
    levels$sample_grouping[levels$sample_grouping=="DKO_12_DKO_1"] <- 12
    levels$sample_grouping[levels$sample_grouping=="DKO_13_DKO_1"] <- 13
    levels$sample_grouping[levels$sample_grouping=="DKO_14_DKO_1"] <- 14
    levels$sample_grouping[levels$sample_grouping=="DKO_15_DKO_4"] <- 15
    levels$sample_grouping[levels$sample_grouping=="DKO_16_DKO_4"] <- 16
    levels$sample_grouping[levels$sample_grouping=="DKO_17_DKO_4"] <- 17
    levels$sample_grouping[levels$sample_grouping=="DKO_18_DKO_4"] <- 18
    levels$sample_grouping[levels$sample_grouping=="DKO_19_DKO_4"] <- 19
    levels$sample_grouping[levels$sample_grouping=="DKO_20_DKO_4"] <- 20
    levels$sample_grouping[levels$sample_grouping=="DKO_21_DKO_4"] <- 21
    levels$sample_grouping[levels$sample_grouping=="DKO_22_DKO_4"] <- 22
    levels$sample_grouping[levels$sample_grouping=="DKO_23_DKO_4"] <- 23
    levels$sample_grouping[levels$sample_grouping=="DKO_24_DKO_5"] <- 24
    levels$sample_grouping[levels$sample_grouping=="DKO_25_DKO_5"] <- 25
    levels$sample_grouping[levels$sample_grouping=="DKO_26_DKO_5"] <- 26
    levels$sample_grouping[levels$sample_grouping=="DKO_27_DKO_3"] <- 27
    levels$sample_grouping[levels$sample_grouping=="DKO_28_DKO_4"] <- 28
    levels$sample_grouping[levels$sample_grouping=="DKO_29_DKO_4"] <- 29
    levels$sample_grouping[levels$sample_grouping=="DKO_30_DKO_4"] <- 30
    levels$sample_grouping[levels$sample_grouping=="DKO_31_DKO_4"] <- 31
    levels$sample_grouping[levels$sample_grouping=="DKO_32_DKO_4"] <- 32
    levels$sample_grouping[levels$sample_grouping=="DKO_33_DKO_6"] <- 33
    levels$sample_grouping[levels$sample_grouping=="DKO_34_DKO_6"] <- 34
    levels$sample_grouping[levels$sample_grouping=="DKO_35_DKO_6"] <- 35
    levels$sample_grouping[levels$sample_grouping=="DKO_36_DKO_6"] <- 36
    levels$sample_grouping[levels$sample_grouping=="DKO_37_DKO_2"] <- 37
    levels$sample_grouping[levels$sample_grouping=="DKO_38_DKO_2"] <- 38
    levels$sample_grouping[levels$sample_grouping=="DKO_39_DKO_2"] <- 39
    levels$sample_grouping[levels$sample_grouping=="DKO_40_DKO_5"] <- 40
    levels$sample_grouping[levels$sample_grouping=="DKO_41_DKO_3"] <- 41
    levels$sample_grouping[levels$sample_grouping=="DKO_42_DKO_3"] <- 42
    levels$sample_grouping[levels$sample_grouping=="DKO_43_DKO_5"] <- 43
    levels$sample_grouping[levels$sample_grouping=="DKO_44_DKO_4"] <- 44
    levels$sample_grouping[levels$sample_grouping=="DKO_49_DKO_2"] <- 45
    levels$sample_grouping[levels$sample_grouping=="DKO_50_DKO_6"] <- 46
    levels$sample_grouping[levels$sample_grouping=="DK_131_DK"] <- 47
    levels$sample_grouping[levels$sample_grouping=="DK_132_DK"] <- 48
    levels$sample_grouping[levels$sample_grouping=="DK_133_DK"] <- 49
    levels$sample_grouping[levels$sample_grouping=="DK_134_DK"] <- 50
    levels$sample_grouping[levels$sample_grouping=="DK_135_DK"] <- 51
    levels$sample_grouping[levels$sample_grouping=="DK_136_DK"] <- 52
    levels$sample_grouping[levels$sample_grouping=="DK_137_DK"] <- 53
    levels$sample_grouping[levels$sample_grouping=="DK_139_DK"] <- 54
    levels$sample_grouping[levels$sample_grouping=="DK_140_DK"] <- 55
    levels$sample_grouping[levels$sample_grouping=="DK_141_DK"] <- 56
    levels$sample_grouping[levels$sample_grouping=="DK_142_DK"] <- 57
    levels$sample_grouping[levels$sample_grouping=="DK_143_DK"] <- 58
    levels$sample_grouping[levels$sample_grouping=="DK_144_DK"] <- 59
    levels$sample_grouping[levels$sample_grouping=="DK_145_DK"] <- 60
    levels$sample_grouping[levels$sample_grouping=="DK_146_DK"] <- 61
    levels$sample_grouping[levels$sample_grouping=="DK_147_DK"] <- 62
    levels$sample_grouping[levels$sample_grouping=="DK_148_DK"] <- 63
    levels$sample_grouping[levels$sample_grouping=="DK_149_DK"] <- 64
    levels$sample_grouping[levels$sample_grouping=="DK_150_DK"] <- 65
    levels$sample_grouping[levels$sample_grouping=="DK_151_DK"] <- 66
    levels$sample_grouping[levels$sample_grouping=="DK_152_DK"] <- 67
    levels$sample_grouping[levels$sample_grouping=="DK_153_DK"] <- 68
    levels$sample_grouping[levels$sample_grouping=="DK_154_DK"] <- 69
    levels$sample_grouping[levels$sample_grouping=="DK_155_DK"] <- 70
    levels$sample_grouping[levels$sample_grouping=="DK_156_DK"] <- 71
    levels$sample_grouping[levels$sample_grouping=="DK_157_DK"] <- 72
    levels$sample_grouping[levels$sample_grouping=="DK_158_DK"] <- 73
    levels$sample_grouping[levels$sample_grouping=="DK_159_DK"] <- 74
    levels$sample_grouping[levels$sample_grouping=="DK_160_DK"] <- 75
    levels$sample_grouping[levels$sample_grouping=="DK_161_DK"] <- 76
    levels$sample_grouping[levels$sample_grouping=="DK_162_DK"] <- 77
    levels$sample_grouping[levels$sample_grouping=="DK_163_DK"] <- 78
    levels$sample_grouping[levels$sample_grouping=="DK_164_DK"] <- 79
    levels$sample_grouping[levels$sample_grouping=="DK_165_DK"] <- 80
    levels$sample_grouping[levels$sample_grouping=="DK_166_DK"] <- 81
    levels$sample_grouping[levels$sample_grouping=="DK_167_DK"] <- 82
    levels$sample_grouping[levels$sample_grouping=="DK_169_DK"] <- 83
    levels$sample_grouping[levels$sample_grouping=="DK_170_DK"] <- 84
    levels$sample_grouping[levels$sample_grouping=="F_51_F"] <- 85
    levels$sample_grouping[levels$sample_grouping=="F_52_F"] <- 86
    levels$sample_grouping[levels$sample_grouping=="F_53_F"] <- 87
    levels$sample_grouping[levels$sample_grouping=="F_54_F"] <- 88
    levels$sample_grouping[levels$sample_grouping=="F_55_F"] <- 89
    levels$sample_grouping[levels$sample_grouping=="F_56_F"] <- 90
    levels$sample_grouping[levels$sample_grouping=="F_57_F"] <- 91
    levels$sample_grouping[levels$sample_grouping=="F_58_F"] <- 92
    levels$sample_grouping[levels$sample_grouping=="F_59_F"] <- 93
    levels$sample_grouping[levels$sample_grouping=="F_60_F"] <- 94
    levels$sample_grouping[levels$sample_grouping=="F_61_F"] <- 95
    levels$sample_grouping[levels$sample_grouping=="F_62_F"] <- 96
    levels$sample_grouping[levels$sample_grouping=="F_63_F"] <- 97
    levels$sample_grouping[levels$sample_grouping=="F_64_F"] <- 98
    levels$sample_grouping[levels$sample_grouping=="F_65_F"] <- 99
    levels$sample_grouping[levels$sample_grouping=="F_66_F"] <- 100
    levels$sample_grouping[levels$sample_grouping=="F_67_F"] <- 101
    levels$sample_grouping[levels$sample_grouping=="F_68_F"] <- 102
    levels$sample_grouping[levels$sample_grouping=="F_69_F"] <- 103
    levels$sample_grouping[levels$sample_grouping=="F_70_F"] <- 104
    levels$sample_grouping[levels$sample_grouping=="F_71_F"] <- 105
    levels$sample_grouping[levels$sample_grouping=="F_72_F"] <- 106
    levels$sample_grouping[levels$sample_grouping=="F_73_F"] <- 107
    levels$sample_grouping[levels$sample_grouping=="F_74_F"] <- 108
    levels$sample_grouping[levels$sample_grouping=="F_79_F"] <- 109
    levels$sample_grouping[levels$sample_grouping=="F_80_F"] <- 110
    levels$sample_grouping[levels$sample_grouping=="F_81_F"] <- 111
    levels$sample_grouping[levels$sample_grouping=="F_82_F"] <- 112
    levels$sample_grouping[levels$sample_grouping=="F_85_F"] <- 113
    levels$sample_grouping[levels$sample_grouping=="F_86_F"] <- 114
    levels$sample_grouping[levels$sample_grouping=="F_87_F"] <- 115
    levels$sample_grouping[levels$sample_grouping=="F_88_F"] <- 116
    levels$sample_grouping[levels$sample_grouping=="F_90_F"] <- 117
    levels$sample_grouping[levels$sample_grouping=="UK_100_UK"] <- 118
    levels$sample_grouping[levels$sample_grouping=="UK_101_UK"] <- 119
    levels$sample_grouping[levels$sample_grouping=="UK_102_UK"] <- 120
    levels$sample_grouping[levels$sample_grouping=="UK_103_UK"] <- 121
    levels$sample_grouping[levels$sample_grouping=="UK_104_UK"] <- 122
    levels$sample_grouping[levels$sample_grouping=="UK_105_UK"] <- 123
    levels$sample_grouping[levels$sample_grouping=="UK_106_UK"] <- 124
    levels$sample_grouping[levels$sample_grouping=="UK_107_UK"] <- 125
    levels$sample_grouping[levels$sample_grouping=="UK_108_UK"] <- 126
    levels$sample_grouping[levels$sample_grouping=="UK_109_UK"] <- 127
    levels$sample_grouping[levels$sample_grouping=="UK_110_UK"] <- 128
    levels$sample_grouping[levels$sample_grouping=="UK_111_UK"] <- 129
    levels$sample_grouping[levels$sample_grouping=="UK_112_UK"] <- 130
    levels$sample_grouping[levels$sample_grouping=="UK_113_UK"] <- 131
    levels$sample_grouping[levels$sample_grouping=="UK_114_UK"] <- 132
    levels$sample_grouping[levels$sample_grouping=="UK_115_UK"] <- 133
    levels$sample_grouping[levels$sample_grouping=="UK_116_UK"] <- 134
    levels$sample_grouping[levels$sample_grouping=="UK_117_UK"] <- 135
    levels$sample_grouping[levels$sample_grouping=="UK_118_UK"] <- 136
    levels$sample_grouping[levels$sample_grouping=="UK_119_UK"] <- 137
    levels$sample_grouping[levels$sample_grouping=="UK_120_UK"] <- 138
    levels$sample_grouping[levels$sample_grouping=="UK_121_UK"] <- 139
    levels$sample_grouping[levels$sample_grouping=="UK_122_UK"] <- 140
    levels$sample_grouping[levels$sample_grouping=="UK_125_UK"] <- 141
    levels$sample_grouping[levels$sample_grouping=="UK_126_UK"] <- 142
    levels$sample_grouping[levels$sample_grouping=="UK_127_UK"] <- 143
    levels$sample_grouping[levels$sample_grouping=="UK_128_UK"] <- 144
    levels$sample_grouping[levels$sample_grouping=="UK_129_UK"] <- 145
    levels$sample_grouping[levels$sample_grouping=="UK_130_UK"] <- 146
    levels$sample_grouping[levels$sample_grouping=="UK_91_UK"] <- 147
    levels$sample_grouping[levels$sample_grouping=="UK_92_UK"] <- 148
    levels$sample_grouping[levels$sample_grouping=="UK_93_UK"] <- 149
    levels$sample_grouping[levels$sample_grouping=="UK_94_UK"] <- 150
    levels$sample_grouping[levels$sample_grouping=="UK_95_UK"] <- 151
    levels$sample_grouping[levels$sample_grouping=="UK_96_UK"] <- 152
    levels$sample_grouping[levels$sample_grouping=="UK_97_UK"] <- 153
    levels$sample_grouping[levels$sample_grouping=="UK_99_UK"] <- 154
    levels$sample_grouping[levels$sample_grouping=="UK_98_UK"] <- 155
    }
    
    levels$sample_grouping <- as.numeric(levels$sample_grouping)
    
  }
  
  # ========================================== #
  #                All samples                 #
  # ========================================== #
  {
    #Calculate Fst values 
    levels2 <- data.frame(levels2[, c(1:3)])
    varcomp.glob(levels,loci, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of management/total (DKOvDLF)
    test.between(loci,rand.unit=levels$grouping,test=levels$sample,nperm=1000) # p-value = 0.727
    # Test for the effect of grouping/total 
    test.between(loci,rand.unit=levels$sample_grouping,test=levels$grouping,nperm=1000) # p-value =  0.001
    # Test for the effect of field_plot/total 
    levels$individual_sample <- row.names(levels)
    test.between(loci,rand.unit=levels$individual_sample,test=levels$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of groupings/management, issue the command
    test.within(loci, test=levels$grouping, within = levels$sample, nperm=1000) # p-value = 0.001
    # To test for the effect of field_plot/management, issue the command
    test.within(loci, test=levels$sample_grouping, within = levels$sample, nperm=1000) # p-value = 0.001
    # To test for the effect of field_plot/grouping, issue the command
    test.within(loci, test=levels$sample_grouping, within = levels$grouping, nperm=1000) # p-value = 0.001
    
  }
  
  # ========================================== #
  #               Field trials                 #
  # ========================================== #
  {
    #Choose field trial samples
    y = nrow(df)
    loci2 <- loci[3700:y,]
    y = nrow(levels)
    levels2 <- levels[3700:y,]
    
    levels2 <- data.frame(levels2[, c(2:3)])
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$sample_grouping, within = levels2$grouping, nperm=1000) # p-value = 0.001
  }
  # ========================================== #
  #                  Organic                   #
  # ========================================== #
  {
    #Choose DKO samples
    y = nrow(df)
    loci2 <- loci[-3700:-y,]
    y = nrow(levels)
    levels2 <- levels[-3700:-y,]
    
    levels2 <- data.frame(levels2[, c(2:3)])
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total 
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.226
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.008
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$sample_grouping, within = levels2$grouping, nperm=1000) # p-value = 0.001
    
  }  
  
  
  
  # ========================================== #
  #      Field trials with clover genotype     #
  # ========================================== #
  {
    #Choose field trial samples
    y = nrow(df)
    loci2 <- loci[3700:y,]
    levels2 <- levels[3700:y,]
    
    #Same number for clover genotypes between field trial sites
    levels2 <- data.frame(levels2$grouping, levels2$field_number, levels2$Same_number_for_DLF, levels2$sample_grouping)
    
    #Different number for clover genotypes between field trial sites
    #levels2 <- data.frame(levels2$grouping, levels2$field_number, levels2$different_between_trials, levels2$sample_grouping)
    
    varcomp.glob(levels2,loci2, diploid=FALSE)
    
    #---- TEST BETWEEN ----#
    # Test for the effect of grouping/total
    test.between(loci2,rand.unit=levels2$sample_grouping,test=levels2$grouping,nperm=1000) # p-value = 0.001
    # Test for the effect of field_plot/total 
    test.between(loci2,rand.unit=levels2$individual_sample,test=levels2$sample_grouping,nperm=1000) # p-value = 0.001
    
    #---- TEST WITHIN ----#
    # To test for the effect of field_plot/within_grouping, issue the command
    test.within(loci2, test=levels2$levels2.sample_grouping, within = levels2$levels2.grouping, nperm=1000) # p-value = 0.001
    # To test for the effect of clovergenotype/within_block, issue the command
    test.within(loci2, test=levels2$levels2.Same_number_for_DLF, within = levels2$levels2.field_number, nperm=1000) # p-value = 0.001
  }
} 

