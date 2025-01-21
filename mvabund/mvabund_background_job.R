#Background mvabund job
# Marissa Wright-La Greca (2024-05-22)

#Clear space
# rm(list=ls())

#Update (2024-06-07)
#Jadi (author) suggested package errors may be occcuring and slowing computation time
#for Composition = TRUE arguement
#Suggests a package reset for the following packages: 

remove.packages("mvabund")
remove.packages("RcppGSL")
remove.packages("Rcpp") 

install.packages("Rcpp") 

install.packages("RcppGSL", repos="http://cran.us.r-project.org")

#packageurl <- "https://cran.r-project.org/src/contrib/Archive/RcppGSL/RcppGSL_0.3.11.tar.gz"
#install.packages(packageurl, repos=NULL, type="source") 

# packageurl <- "https://cran.r-project.org/src/contrib/Archive/mvabund/mvabund_3.12.3.tar.gz"
# install.packages(packageurl, repos=NULL, type="source") 

install.packages("mvabund")

#Packages ----
# install.packages(tidyverse)
# install.packages(vegan)
# install.packages(mvabund)
# install.packages("tibble")
# install.packages(RColorBrewer) 
# install.packages(phyloseq)
# install.packages(microbiome)
# install.packages(parallel)
# install.packages("mvabund")
# install.packages("dada2") #changed from data2, assuming typo
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("dada2", force = TRUE)

#Library ----
library(tidyverse)
#library(data2) #changed from dada2
library(mvabund)
#library(RColorBrewer) 
library(phyloseq)
library(microbiome)
library(tibble)
library(ggplot2)
library(reshape2)


#library(parallel) #to run parallel processing for manyglm fxn

#Load phyloseq object ----

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current.path)

# User set variables
# data.FN <- "~/Documents/Marissa_Github/GitHub/mb2021_phyloseq/Old RDS files/mb2021_filtered_NOT_rarefied.rds" # MacPro (CSR)
# #data.FN <- "" # MacPro (CSR)

data.FN <- "~Phyloseq and microbiome analysis/Old RDS files/mb2021_filtered_NOT_rarefied.rds"
data.FN <- "~Phyloseq and microbiome analysis/Old RDS files/MU42022_filtered_NOT_rarefied.rds"

# Load data
# mb2021_filtered_NOT_rarefied <- readRDS(data.FN)

MU42022_filtered_NOT_rarefied <- readRDS(data.FN)

# Rename as per tutorial
pseq <- MU42022_filtered_NOT_rarefied
pseq <- MU42022_filtered_Oct92024
pseq <- MU42022_spat_unrarefied
pseq <- MU42022_normalized_relative

# Subset samples ----
#mb2021 project
pseq <- subset_samples(pseq, Age %in% c("1 dpf"))  
#pseq <- subset_samples(pseq, Age %in% c("Spat"))
#pseq <- subset_samples(pseq, !Family %in% c("9"))
#pseq <- subset_samples(pseq, !Treatment %in% c("Low salinity"))

#MU42022 project
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")
pseq <- subset_samples(pseq, Age %in% c("Spat"))

#Sanity check
#check if any OTUs are not present in any samples (want false)
any(taxa_sums(pseq) == 0)

#if true
pseq_filtered <- prune_taxa(taxa_sums(pseq) > 0, pseq)
any(taxa_sums(pseq_filtered) == 0)

pseq <- pseq_filtered


# Agglomerating family names
# pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)]  <- 1
# pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
# pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
# pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4
# pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

# Subset option for agglomerated family names
#pseq <- subset_samples(pseq, !Family %in% c("1"))


# Subset for specific ASVs ----
# look at certain ASVs found to be significant but reduce treatments 

# Designate day 1 signif ASVs
# testing
taxa_to_subset <- c("ASV21", "ASV34", "ASV40", "ASV42", "ASV44", "ASV48", "ASV51", "ASV60", "ASV76", "ASV79", "ASV90", "ASV91", "ASV92", "ASV100", "ASV109", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV130", "ASV131", "ASV150", "ASV163", "ASV176", "ASV181", "ASV187", "ASV189", "ASV193", "ASV195", "ASV216", "ASV218", "ASV219", "ASV237", "ASV240", "ASV244", "ASV257", "ASV262", "ASV309", "ASV313", "ASV315", "ASV319", "ASV326", "ASV333", "ASV337", "ASV343", "ASV349", "ASV351", "ASV356", "ASV375", "ASV379", "ASV385", "ASV387", "ASV394", "ASV397", "ASV401", "ASV416", "ASV417", "ASV428", "ASV431", "ASV447", "ASV450", "ASV468", "ASV484", "ASV487", "ASV499", "ASV501", "ASV510", "ASV523", "ASV524", "ASV538", "ASV556", "ASV559", "ASV563", "ASV575", "ASV579", "ASV595", "ASV608", "ASV620", "ASV641", "ASV649", "ASV652", "ASV703", "ASV738", "ASV756", "ASV758", "ASV789", "ASV798", "ASV818", "ASV819", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV874", "ASV882", "ASV889", "ASV916", "ASV918", "ASV919", "ASV929", "ASV952", "ASV955", "ASV960", "ASV982", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1030", "ASV1035", "ASV1036", "ASV1061", "ASV1070", "ASV1074", "ASV1078", "ASV1082", "ASV1086", "ASV1110", "ASV1111", "ASV1138", "ASV1146", "ASV1153", "ASV1170")
#taxa_to_subset <- taxa_to_subset[1:11]
#taxa_to_subset <- taxa_to_subset[1:50]
taxa_to_subset <- taxa_to_subset[1:100]
#taxa_to_subset <- taxa_to_subset[1:200] #NOT ENOUGH TAXA
#taxa_to_subset <- taxa_to_subset[1:400] #NOT ENOUGH TAXA


# Spat significant ASVs (without a p-value adjustment)
# taxa_to_subset <- c("ASV17",  "ASV30",   "ASV37",   "ASV49",   "ASV216",  "ASV237",  "ASV240",  "ASV254",  "ASV271",  "ASV314", "ASV325", "ASV347"
#                              , "ASV348", "ASV383", "ASV391", "ASV412", "ASV446", "ASV507", "ASV553", "ASV571", "ASV580", "ASV595", "ASV656", "ASV751"
#                              , "ASV841",  "ASV845",  "ASV900",  "ASV927",  "ASV964",  "ASV1088", "ASV1162", "ASV1164"
#                             )

# 1 dpf signif ASVs
#taxa_to_subset <- c("ASV21", "ASV34", "ASV40", "ASV42", "ASV44", "ASV48", "ASV51", "ASV60", "ASV76", "ASV79", "ASV90", "ASV91", "ASV92", "ASV100", "ASV109", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV130", "ASV131", "ASV150", "ASV163", "ASV176", "ASV181", "ASV187", "ASV189", "ASV193", "ASV195", "ASV216", "ASV218", "ASV219", "ASV237", "ASV240", "ASV244", "ASV257", "ASV262", "ASV309", "ASV313", "ASV315", "ASV319", "ASV326", "ASV333", "ASV337", "ASV343", "ASV349", "ASV351", "ASV356", "ASV375", "ASV379", "ASV385", "ASV387", "ASV394", "ASV397", "ASV401", "ASV416", "ASV417", "ASV428", "ASV431", "ASV447", "ASV450", "ASV468", "ASV484", "ASV487", "ASV499", "ASV501", "ASV510", "ASV523", "ASV524", "ASV538", "ASV556", "ASV559", "ASV563", "ASV575", "ASV579", "ASV595", "ASV608", "ASV620", "ASV641", "ASV649", "ASV652", "ASV703", "ASV738", "ASV756", "ASV758", "ASV789", "ASV798", "ASV818", "ASV819", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV874", "ASV882", "ASV889", "ASV916", "ASV918", "ASV919", "ASV929", "ASV952", "ASV955", "ASV960", "ASV982", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1030", "ASV1035", "ASV1036", "ASV1061", "ASV1070", "ASV1074", "ASV1078", "ASV1082", "ASV1086", "ASV1110", "ASV1111", "ASV1138", "ASV1146", "ASV1153", "ASV1170")


# Prune taxa (optional)
# Reporting
# print(paste0("Pruning to only retain the following taxa/ASVs: "))
# print(taxa_to_subset)
# print(paste0("There are ", length(taxa_to_subset), " taxa in this analysis"))
# 
# # Prune taxa
pseq <- prune_taxa(taxa_to_subset, pseq)

#Use relative abundance? = Does not help with zero-inflation
#pseq <- microbiome::transform(pseq, "compositional")


# Create mvabund object from phyloseq object ----
   # note: not sure why but metadata must first be transformed into matrix then data frame to convert properly
fact1 = sample_data(pseq)
fact = as.matrix.data.frame(fact1)
fact = as.data.frame(fact)
View(fact)

#Make treatments into factors rather than characters
fact$Sample.ID <- as.factor(fact$Sample.ID)
# fact$Library_Name <- as.factor(fact$Library_Name)
fact$Treatment <- as.factor(fact$Treatment)
fact$Genetics <- as.factor(fact$Genetics)
fact$Age <- as.factor(fact$Age)

# Reset row names to be numbered
rownames(fact) <- NULL

ASV_data_cleaned <- pseq@otu_table
ASV_data_cleaned <- as.data.frame(pseq@otu_table)

#Remove ASVs with sum of zero across all samples
ASV_data_cleaned <- ASV_data_cleaned[, colSums(ASV_data_cleaned) > 0]

#Checking for errors in data - no NAs
dim(ASV_data_cleaned)

any(is.na(rownames(ASV_data_cleaned))) #FALSE
anyNA(ASV_data_cleaned) #FALSE

#fact$numberReads <- rowSums(ASV_data_cleaned)
#dat_mvabund <- mvabund(ASV_data_cleaned, row.names=0, check.row = TRUE, check.names = TRUE, var.names = "") #row.names argument inputs row names (sample info) into mvabund object
dat_mvabund <- mvabund(ASV_data_cleaned)
View(dat_mvabund) #Note: Check row names correct (contains sample IDs)

#Check if converted to mvabund object properly
is.mvabund(dat_mvabund)

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

ASV_data_cleaned <- as.data.frame(pseq@otu_table)
rowSums(ASV_data_cleaned)
fact$numberReads <- rowSums(ASV_data_cleaned)
dat_mvabund <- mvabund(ASV_data_cleaned)

# Standard method (log number of reads)
Mod_log <- manyglm(dat_mvabund ~ Treatment * Genetics + offset(log(numberReads)), family = "negative.binomial", data = fact)
plot(Mod_log)
summary(Mod_log)
#State time
Sys.time()

# Standard method, setting composition argument
#Note: Comp = TRUE may not allow for p adjustment (see pg 12: https://cran.r-project.org/web/packages/mvabund/mvabund.pdf)
dat_nb_compositionT_Spat_LS = manyglm(dat_mvabund ~ Treatment, family="negative.binomial", data = fact, composition=TRUE, show.warning = TRUE)
dat_nb_Spat_LS = manyglm(dat_mvabund ~ Treatment * Genetics, family="negative.binomial", data = fact, composition=FALSE)
plot(dat_nb_compositionT_Spat_LS) #check residuals to see model fit

summary(dat_nb_Spat_LS)
plot(dat_nb_Spat_LS)

# Standardized residuals
std_resid <- residuals(dat_nb_compositionT_Spat_LS, type = "standardized")
hist(std_resid, main = "Histogram of Standardized Residuals")

# Checking for patterns in residuals
plot(fitted(dat_nb_compositionT_Spat_LS), std_resid, 
     xlab = "Fitted values", ylab = "Standardized Residuals")
abline(h = 0, col = "red")


#Visualize data some more

#Histograms for subset of ASVs
# Step 1: Create a 'SampleID' column from the row names
ASV_subset$SampleID <- rownames(ASV_subset)

# Step 2: Melt the data into long format
ASV_long <- melt(ASV_subset, id.vars = "SampleID")

# Step 3: Check the first few rows to confirm the transformation
head(ASV_long)

# Create histograms for each ASV
ggplot(ASV_long, aes(x = value)) +
  geom_histogram(binwidth = 10, fill = "skyblue", color = "black") +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Distribution of ASV Abundances", x = "Abundance", y = "Frequency") +
  theme_minimal()

# Count number of zeros per ASV
zero_counts <- colSums(ASV_data_cleaned == 0)

# Bar plot of zero counts per ASV
zero_df <- data.frame(ASV = names(zero_counts), ZeroCount = zero_counts)
ggplot(zero_df, aes(x = ASV, y = ZeroCount)) +
  geom_bar(stat = "identity", fill = "lightcoral") +
  labs(title = "Number of Zeros in Each ASV", x = "ASV", y = "Zero Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Density plots
ASV_long_nonzero <- ASV_long[ASV_long$value > 0, ]
ggplot(ASV_long_nonzero, aes(x = value)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Density Plot of ASV Abundances (Non-Zero Values)", x = "Abundance", y = "Density") +
  theme_minimal()





Sys.time()
dat_aov_unadj_Spat_LS <- anova(dat_nb_compositionT_Spat_LS, p.uni = "unadjusted",show.time = "all", show.warning = TRUE)
Sys.time()

#composition = TRUE (round 2)
#dat_nb_compositionT_Spat_LS_round2 = manyglm(dat_mvabund ~ Treatment * Family, family="negative.binomial", data = fact, composition=TRUE)
dat_nb_compositionT_Spat_LS_round2 = manyglm(dat_mvabund ~ Treatment * Family, family="negative.binomial", data = fact, composition=TRUE)
Sys.time()

dat_aov_unadj_Spat_LS_round2 <- anova(dat_nb_compositionT_Spat_LS_round2, p.uni = "unadjusted", show.time = "all")
Sys.time()

#anova is taking a long time - model not converging properly?

#Try a zero-inflated model: To handle zeros which may be causing model to not converge
install.packages("glmmTMB")
library(glmmTMB)
zinb_model <- glmmTMB(dat_mvabund ~ Treatment * Genetics, 
                      data = fact, 
                      family = nbinom2, 
                      ziformula = ~ 1) 





# Run ANOVA ----
# note: need to select adjusted or unadjusted p
dat_aov_unadj_Spat_LS <- anova(dat_nb_compositionT_Spat_LS, p.uni = "unadjusted")
dat_aov_unadj_Spat_LS


### Result exploration #

#Select alpha level----
length(which(dat_aov_unadj_Spat_LS_round2$uni.p[2,] < 0.05)) # Number of ASVs affected by treatment
length(which(dat_aov_unadj_Spat_LS_round2$uni.p[3,] < 0.05)) # Number of ASVs affected by family
length(which(dat_aov_unadj_Spat_LS_round2$uni.p[4,] < 0.05)) # Number of ASVs affected by interaction

#Select ASVs that are different between factors ----
diff_05_treat_unadj_Spat <- which(dat_aov_unadj_Spat_LS_round2$uni.p[2,] < 0.05) # [2] is referring to the first factor of my model , in my case "Treatment"
diff_05_family_unadj_Spat <- which(dat_aov_unadj_Spat_LS_round2$uni.p[3,] < 0.05) # [3] is referring to the second factor of my model , in my case "Family" ***day time but thats wrong
diff_05_interact_unadj_Spat <- which(dat_aov_unadj_Spat_LS_round2$uni.p[4,] < 0.05) # [4] is referring to the third factor of my model , in my case "interaction Treatment and Family"

#ASV names ----
names_diff_treat_unadj_Spat <- names(diff_05_treat_unadj_Spat)
names_diff_family_unadj_Spat <- names(diff_05_family_unadj_Spat)
names_diff_interact_unadj_Spat <- names(diff_05_interact_unadj_Spat)

#Store results ----
# TODO





