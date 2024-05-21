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
#pseq <- subset_samples(pseq, Age %in% c("1 dpf"))
pseq <- subset_samples(pseq, Age %in% c("Spat"))
pseq <- subset_samples(pseq, !Family %in% c("9"))
pseq <- subset_samples(pseq, !Treatment %in% c("High salinity"))

#correct family names
pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

#pseq <- subset_samples(pseq, !Family %in% c("1"))

#Look at certain ASVs found to be significant but reduce treatments 

#Spat signif ASVs
Signif_ASVs_Spat_unadj <- c("ASV17",  "ASV30",   "ASV37",   "ASV49",   "ASV216",  "ASV237",  "ASV240",  "ASV254",  "ASV271",  "ASV314", "ASV325", "ASV347", 
                    "ASV348", "ASV383", "ASV391", "ASV412", "ASV446", "ASV507", "ASV553", "ASV571", "ASV580", "ASV595", "ASV656", "ASV751", 
                    "ASV841",  "ASV845",  "ASV900",  "ASV927",  "ASV964",  "ASV1088", "ASV1162", "ASV1164")

#1 dpf signif ASVs
#Signif_ASVs_1dpf_unadj <- c("ASV21", "ASV34", "ASV40", "ASV42", "ASV44", "ASV48", "ASV51", "ASV60", "ASV76", "ASV79", "ASV90", "ASV91", "ASV92", "ASV100", "ASV109", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV130", "ASV131", "ASV150", "ASV163", "ASV176", "ASV181", "ASV187", "ASV189", "ASV193", "ASV195", "ASV216", "ASV218", "ASV219", "ASV237", "ASV240", "ASV244", "ASV257", "ASV262", "ASV309", "ASV313", "ASV315", "ASV319", "ASV326", "ASV333", "ASV337", "ASV343", "ASV349", "ASV351", "ASV356", "ASV375", "ASV379", "ASV385", "ASV387", "ASV394", "ASV397", "ASV401", "ASV416", "ASV417", "ASV428", "ASV431", "ASV447", "ASV450", "ASV468", "ASV484", "ASV487", "ASV499", "ASV501", "ASV510", "ASV523", "ASV524", "ASV538", "ASV556", "ASV559", "ASV563", "ASV575", "ASV579", "ASV595", "ASV608", "ASV620", "ASV641", "ASV649", "ASV652", "ASV703", "ASV738", "ASV756", "ASV758", "ASV789", "ASV798", "ASV818", "ASV819", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV874", "ASV882", "ASV889", "ASV916", "ASV918", "ASV919", "ASV929", "ASV952", "ASV955", "ASV960", "ASV982", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1030", "ASV1035", "ASV1036", "ASV1061", "ASV1070", "ASV1074", "ASV1078", "ASV1082", "ASV1086", "ASV1110", "ASV1111", "ASV1138", "ASV1146", "ASV1153", "ASV1170")


pseq <- prune_taxa(Signif_ASVs_Spat_unadj, pseq)

fact1 = sample_data(pseq)
fact = as.matrix.data.frame(fact1)
fact = as.data.frame(fact)
ASV_data_cleaned <- pseq@otu_table

ASV_data_cleaned <- as.data.frame(pseq@otu_table)
rowSums(ASV_data_cleaned)
fact$numberReads <- rowSums(ASV_data_cleaned)
dat_mvabund <- mvabund(ASV_data_cleaned)

#dat_nb_Spat_HS <- manyglm(dat_mvabund ~ Treatment * Family + offset(log(numberReads)), family = "negative.binomial", data = fact)
dat_nb_compositionT_Spat_LS = manyglm(dat_mvabund ~ Treatment * Family, family="negative.binomial", data = fact, composition=TRUE)


dat_aov_unadj_Spat_LS <- anova(dat_nb_compositionT_Spat_LS, p.uni = "unadjusted") 

length(which(dat_aov_unadj_Spat_LS$uni.p[2,] < 0.05)) # Number of ASVs affected by treatment
length(which(dat_aov_unadj_Spat_LS$uni.p[3,] < 0.05)) # Number of ASVs affected by family
length(which(dat_aov_unadj_Spat_LS$uni.p[4,] < 0.05)) # Number of ASVs affected by interaction

diff_05_treat_unadj_Spat_LS <- which(dat_aov_unadj_Spat_LS$uni.p[2,] < 0.05) # [2] is referring to the first factor of my model , in my case "Treatment"
diff_05_family_unadj_Spat_LS <- which(dat_aov_unadj_Spat_LS$uni.p[3,] < 0.05) # [3] is referring to the second factor of my model , in my case "Family" ***day time but thats wrong
diff_05_interact_unadj_Spat_LS <- which(dat_aov_unadj_Spat_LS$uni.p[4,] < 0.05) # [4] is referring to the third factor of my model , in my case "interaction Treatment and Family"

names_diff_treat_unadj_Spat_LS <- names(diff_05_treat_unadj_Spat_LS)
names_diff_family_unadj_Spat_LS <- names(diff_05_family_unadj_Spat_LS)
names_diff_interact_unadj_Spat_LS <- names(diff_05_interact_unadj_Spat_LS)







