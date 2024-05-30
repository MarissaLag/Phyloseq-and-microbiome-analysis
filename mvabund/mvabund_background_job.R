#Background mvabund job
# Marissa Wright-La Greca (2024-05-22)

#Clear space
# rm(list=ls())

#Library ----
library(tidyverse)
library(vegan)
library(dada2) 
library(mvabund)
library(RColorBrewer) 
library(phyloseq)
library(microbiome)
library(parallel) #to run parallel processing for manyglm fxn

#Load phyloseq object ----

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current.path)

# User set variables
data.FN <- "~/Documents/Marissa_Github/GitHub/mb2021_phyloseq/Old RDS files/mb2021_filtered_NOT_rarefied.rds" # MacPro (CSR)
#data.FN <- "" # MacPro (CSR)

# Load data
mb2021_filtered_NOT_rarefied <- readRDS(data.FN)

# Rename as per tutorial
pseq <- mb2021_filtered_NOT_rarefied

# Subset samples ----
pseq <- subset_samples(pseq, Age %in% c("1 dpf"))   # retain only 1 DPF samples
#pseq <- subset_samples(pseq, Age %in% c("Spat"))
#pseq <- subset_samples(pseq, !Family %in% c("9"))
#pseq <- subset_samples(pseq, !Treatment %in% c("Low salinity"))

# Agglomerating family names
pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)]  <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4
pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

# Subset option for agglomerated family names
#pseq <- subset_samples(pseq, !Family %in% c("1"))


# Subset for specific ASVs ----
# look at certain ASVs found to be significant but reduce treatments 

# Designate spat signif ASVs
# testing
taxa_to_subset <- c("ASV21", "ASV34", "ASV40", "ASV42", "ASV44", "ASV48", "ASV51", "ASV60", "ASV76", "ASV79", "ASV90", "ASV91", "ASV92", "ASV100", "ASV109", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV130", "ASV131", "ASV150", "ASV163", "ASV176", "ASV181", "ASV187", "ASV189", "ASV193", "ASV195", "ASV216", "ASV218", "ASV219", "ASV237", "ASV240", "ASV244", "ASV257", "ASV262", "ASV309", "ASV313", "ASV315", "ASV319", "ASV326", "ASV333", "ASV337", "ASV343", "ASV349", "ASV351", "ASV356", "ASV375", "ASV379", "ASV385", "ASV387", "ASV394", "ASV397", "ASV401", "ASV416", "ASV417", "ASV428", "ASV431", "ASV447", "ASV450", "ASV468", "ASV484", "ASV487", "ASV499", "ASV501", "ASV510", "ASV523", "ASV524", "ASV538", "ASV556", "ASV559", "ASV563", "ASV575", "ASV579", "ASV595", "ASV608", "ASV620", "ASV641", "ASV649", "ASV652", "ASV703", "ASV738", "ASV756", "ASV758", "ASV789", "ASV798", "ASV818", "ASV819", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV874", "ASV882", "ASV889", "ASV916", "ASV918", "ASV919", "ASV929", "ASV952", "ASV955", "ASV960", "ASV982", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1030", "ASV1035", "ASV1036", "ASV1061", "ASV1070", "ASV1074", "ASV1078", "ASV1082", "ASV1086", "ASV1110", "ASV1111", "ASV1138", "ASV1146", "ASV1153", "ASV1170")
#taxa_to_subset <- taxa_to_subset[1:11]
#taxa_to_subset <- taxa_to_subset[1:50]
#taxa_to_subset <- taxa_to_subset[1:100]
#taxa_to_subset <- taxa_to_subset[1:200]
taxa_to_subset <- taxa_to_subset[1:400]


# Spat significant ASVs (without a p-value adjustment)
# taxa_to_subset <- c("ASV17",  "ASV30",   "ASV37",   "ASV49",   "ASV216",  "ASV237",  "ASV240",  "ASV254",  "ASV271",  "ASV314", "ASV325", "ASV347"
#                              , "ASV348", "ASV383", "ASV391", "ASV412", "ASV446", "ASV507", "ASV553", "ASV571", "ASV580", "ASV595", "ASV656", "ASV751"
#                              , "ASV841",  "ASV845",  "ASV900",  "ASV927",  "ASV964",  "ASV1088", "ASV1162", "ASV1164"
#                             )

# 1 dpf signif ASVs
#taxa_to_subset <- c("ASV21", "ASV34", "ASV40", "ASV42", "ASV44", "ASV48", "ASV51", "ASV60", "ASV76", "ASV79", "ASV90", "ASV91", "ASV92", "ASV100", "ASV109", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV130", "ASV131", "ASV150", "ASV163", "ASV176", "ASV181", "ASV187", "ASV189", "ASV193", "ASV195", "ASV216", "ASV218", "ASV219", "ASV237", "ASV240", "ASV244", "ASV257", "ASV262", "ASV309", "ASV313", "ASV315", "ASV319", "ASV326", "ASV333", "ASV337", "ASV343", "ASV349", "ASV351", "ASV356", "ASV375", "ASV379", "ASV385", "ASV387", "ASV394", "ASV397", "ASV401", "ASV416", "ASV417", "ASV428", "ASV431", "ASV447", "ASV450", "ASV468", "ASV484", "ASV487", "ASV499", "ASV501", "ASV510", "ASV523", "ASV524", "ASV538", "ASV556", "ASV559", "ASV563", "ASV575", "ASV579", "ASV595", "ASV608", "ASV620", "ASV641", "ASV649", "ASV652", "ASV703", "ASV738", "ASV756", "ASV758", "ASV789", "ASV798", "ASV818", "ASV819", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV874", "ASV882", "ASV889", "ASV916", "ASV918", "ASV919", "ASV929", "ASV952", "ASV955", "ASV960", "ASV982", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1030", "ASV1035", "ASV1036", "ASV1061", "ASV1070", "ASV1074", "ASV1078", "ASV1082", "ASV1086", "ASV1110", "ASV1111", "ASV1138", "ASV1146", "ASV1153", "ASV1170")


# Prune taxa (optional)
# Reporting
print(paste0("Pruning to only retain the following taxa/ASVs: "))
print(taxa_to_subset)
print(paste0("There are ", length(taxa_to_subset), " taxa in this analysis"))

# Prune taxa
pseq <- prune_taxa(taxa_to_subset, pseq)


# Create mvabund object from phyloseq object ----
   # note: not sure why but metadata must first be transformed into matrix then data frame to convert properly
fact1 = sample_data(pseq)
fact = as.matrix.data.frame(fact1)
fact = as.data.frame(fact)

ASV_data_cleaned <- pseq@otu_table
ASV_data_cleaned <- as.data.frame(pseq@otu_table)
#rowSums(ASV_data_cleaned)

fact$numberReads <- rowSums(ASV_data_cleaned)
dat_mvabund <- mvabund(ASV_data_cleaned)

#Run manyglms ----
# <<<<<<< Updated upstream
# #VERY slow - ways to speed it up? Trying parallel processing
# 
# # Define the number of cores you want to use
# num_cores <- detectCores() - 1  # Use one less than the total number of cores
# 
# # Create a cluster
# cl <- makeCluster(num_cores)
# 
# # Define the model and anova functions
# run_model <- function(dat_mvabund, fact) {
#   manyglm(dat_mvabund ~ Treatment * Family, family = "negative.binomial", data = fact, composition = TRUE)
# }
# 
# run_anova <- function(run_model) {
#   anova(run_model, p.uni = "unadjusted")
# }
# 
# # Export necessary objects and libraries to the cluster
# clusterExport(cl, c("dat_mvabund", "fact", "run_model", "run_anova"))
# clusterEvalQ(cl, library(mvabund))
# 
# # Run the model and ANOVA in parallel
# results <- parLapply(cl, 1:num_cores, function(x) {
#   model <- run_model(dat_mvabund, fact)
#   anova_result <- run_anova(model)
#   return(anova_result)
# # })
# 
# # Stop the cluster
# stopCluster(cl)
# 
# # Combine or analyze the results as needed
# # Here, we're assuming that the first result is representative (if they are all the same)
# dat_aov_unadj_Spat <- results[[1]]
# =======
# # note: VERY slow - ways to speed it up? Trying parallel processing
# >>>>>>> Stashed changes


#### Parallel processing method (in development) ####
# # Define the number of cores you want to use
# num_cores <- detectCores() - 1  # Use one less than the total number of cores
# 
# # Create a cluster
# cl <- makeCluster(num_cores)
# 
# # Define the model and anova functions
# run_model <- function(dat_mvabund, fact) {
#   manyglm(dat_mvabund ~ Treatment * Family, family = "negative.binomial", data = fact, composition = TRUE)
# }
# 
# run_anova <- function(run_model) {
#   anova(run_model, p.uni = "unadjusted")
# }
# 
# #objects and libraries to the cluster
# clusterExport(cl, c("dat_mvabund", "fact", "run_model", "run_anova"))
# clusterEvalQ(cl, library(mvabund))
# 
# # Run the model and anova in parallel
# results <- parLapply(cl, 1:num_cores, function(x) {
#   model <- run_model(dat_mvabund, fact)
#   run_anova(run_model) #error here?
# })
# 
# #Should code be like this
# results <- parLapply(cl, 1:num_cores, function(x) {
#   model <- run_model(dat_mvabund, fact)
#   run_anova(model)
# })
# 
# results <- parLapply(cl, 1:num_cores, function(x) {
#   run_model(dat_mvabund, fact)
# })
# 
# 
# stopCluster(cl)
# 
# dat_aov_unadj_Spat <- results[[1]]
#### /END/ Parallel processing method (in development) ####


# Standard method
#dat_nb_Spat_HS <- manyglm(dat_mvabund ~ Treatment * Family + offset(log(numberReads)), family = "negative.binomial", data = fact)

Sys.time()
# Standard method, setting composition argument
dat_nb_compositionT_Spat_LS = manyglm(dat_mvabund ~ Treatment * Family, family="negative.binomial", data = fact, composition=FALSE)
Sys.time()
dat_aov_unadj_Spat_LS <- anova(dat_nb_compositionT_Spat_LS, p.uni = "unadjusted")
Sys.time()

dat_nb_compositionT_Spat_LS_round2 = manyglm(dat_mvabund ~ Treatment * Family, family="negative.binomial", data = fact, composition=TRUE)
Sys.time()
dat_aov_unadj_Spat_LS_round2 <- anova(dat_nb_compositionT_Spat_LS_round2, p.uni = "unadjusted")
Sys.time()

# Run ANOVA ----
# note: need to select adjusted or unadjusted p
dat_aov_unadj_Spat_LS <- anova(dat_nb_compositionT_Spat_LS, p.uni = "unadjusted")
dat_aov_unadj_Spat_LS


### Result exploration #

#Select alpha level----
length(which(dat_aov_unadj_Spat$uni.p[2,] < 0.05)) # Number of ASVs affected by treatment
length(which(dat_aov_unadj_Spat$uni.p[3,] < 0.05)) # Number of ASVs affected by family
length(which(dat_aov_unadj_Spat$uni.p[4,] < 0.05)) # Number of ASVs affected by interaction

#Select ASVs that are different between factors ----
diff_05_treat_unadj_Spat <- which(dat_aov_unadj_Spat$uni.p[2,] < 0.05) # [2] is referring to the first factor of my model , in my case "Treatment"
diff_05_family_unadj_Spat <- which(dat_aov_unadj_Spat$uni.p[3,] < 0.05) # [3] is referring to the second factor of my model , in my case "Family" ***day time but thats wrong
diff_05_interact_unadj_Spat <- which(dat_aov_unadj_Spat$uni.p[4,] < 0.05) # [4] is referring to the third factor of my model , in my case "interaction Treatment and Family"

#ASV names ----
names_diff_treat_unadj_Spat <- names(diff_05_treat_unadj_Spat)
names_diff_family_unadj_Spat <- names(diff_05_family_unadj_Spat)
names_diff_interact_unadj_Spat <- names(diff_05_interact_unadj_Spat)

#Store results ----
# TODO





