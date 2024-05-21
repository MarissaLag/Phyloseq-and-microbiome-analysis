#Background mvabund job

library(tidyverse)
library(vegan)
library(dada2) #is this supposed to be dada2 (was data2)??
library(mvabund)
library(RColorBrewer) 
library(phyloseq)
library(microbiome)

setwd("/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files")

mb2021_filtered_NOT_rarefied <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/mb2021_filtered_NOT_rarefied.rds")

pseq <- mb2021_filtered_NOT_rarefied

pseq <- subset_samples(pseq, Age %in% c("Spat"))

#correct family names
pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

pseq <- subset_samples(pseq, !Family %in% c("1"))

fact1 = sample_data(pseq)
fact = as.matrix.data.frame(fact1)
fact = as.data.frame(fact)
ASV_data_cleaned <- pseq@otu_table

ASV_data_cleaned <- as.data.frame(pseq@otu_table)
rowSums(ASV_data_cleaned)
fact$numberReads <- rowSums(ASV_data_cleaned)
dat_mvabund <- mvabund(ASV_data_cleaned)

dat_nb_spat <- manyglm(dat_mvabund ~ Treatment * Family + offset(log(numberReads)), family = "negative.binomial", data = fact)

dat_aov_unadj_spat_F1_removed <- anova(dat_nb_spat, p.uni = "unadjusted") 

length(which(dat_aov_unadj_spat_F1_removed$uni.p[2,] < 0.05)) # Number of ASVs affected by treatment
length(which(dat_aov_unadj_spat_F1_removed$uni.p[3,] < 0.05)) # Number of ASVs affected by family
length(which(dat_aov_unadj_spat_F1_removed$uni.p[4,] < 0.05)) # Number of ASVs affected by interaction

diff_05_treat_unadj_spat_F1_removed <- which(dat_aov_unadj_spat_F1_removed$uni.p[2,] < 0.05) # [2] is referring to the first factor of my model , in my case "Treatment"
diff_05_time_unadj_spat_F1_removed <- which(dat_aov_unadj_spat_F1_removed$uni.p[3,] < 0.05) # [3] is referring to the second factor of my model , in my case "Family" ***day time but thats wrong
diff_05_interact_unadj_spat_F1_removed <- which(dat_aov_unadj_spat_F1_removed$uni.p[4,] < 0.05) # [4] is referring to the third factor of my model , in my case "interaction Treatment and Family"

names_diff_treat_unadj_spat_F1_removed <- names(diff_05_treat_unadj_spat_F1_removed)
names_diff_family_unadj_spat_F1_removed <- names(diff_05_time_unadj_spat_F1_removed)
names_diff_interact_unadj_spat_F1_removed <- names(diff_05_interact_unadj_spat_F1_removed)



