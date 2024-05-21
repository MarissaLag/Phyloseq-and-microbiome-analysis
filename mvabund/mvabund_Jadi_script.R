#Using mvabund to run multivariate analysis

#mvabund uses GLMs to model abundance variation for each species in the data set
#This makes the analysis more powerful than other multivariate tests such as PERMANOVA that
#does nothing to account for high variablility in the data set

# More info on mvabund here:
# https://environmentalcomputing.net/statistics/mvabund/
#MVABUND involves very long steps and different OPTIONS for the RUN. Do not run the entire script. You need to chose the best options for your study first.

install.packages("mvabund")
install.packages("data2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2", force = TRUE)

pseq <- Marissa_mb2021_filtered_20240203

#trying with unrarefied data too (has been filtered already)

pseq <- mb2021_filtered_NOT_rarefied

pseq <- subset_samples(pseq, Age %in% c("Spat"))

#correct family names
pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

library(tidyverse)
library(vegan)
library(dada2) #is this supposed to be dada2 (was data2)??
library(mvabund)
library(RColorBrewer) 


# Load your factor table , in my case is "fact"
#For some reason with sam_data in phyloseq - to convert to a data fram you must 1st convert to a matrix (may be a software bug?)
fact1 = sam_data(pseq)
fact = as.matrix.data.frame(fact1)
fact = as.data.frame(fact)
str(fact)
View(fact)
head(fact)

# Load your data table (with counts) in my case is "ASV_data_cleaned"
# Load data (samples in rows and ASVs in column) and metadata

ASV_data_cleaned <- pseq@otu_table

ASV_data_cleaned <- as.data.frame(pseq@otu_table)
str(ASV_data_cleaned)
# #### Multivariate analyses #### In my case, I wanted to test the effect of Treatment and Time. You will have to replace these two with Sample_type or Host_species, depends on what you want to test. 

#PERMANOVA and PERMDISP tests are here just to compare power between these tests and mvabund 
# # PERMANOVA
adonis2(ASV_data_cleaned ~ fact$Treatment*fact$Age*fact$Family)

#It is also worth running PERMDISP
# In fact, I would always run this alongside PERMANOVA to determine if your effects are on the location, dispersion or both (see Warton et al. 2012, Methods in Ecol. & Evol.)

dat_disp = betadisper(dat_bc, fact$treatment, type = "centroid")
par(pty = 's')
boxplot(dat_disp)

permutest(dat_disp) #

dat_disp = betadisper(dat_bc, fact$time, type = "centroid")
par(pty = 's')
boxplot(dat_disp)

permutest(dat_disp) # 


#****************multivariate GLM******************
dat_mvabund <- mvabund(ASV_data_cleaned)

# # GLM multivariate model
# # incorporate total number of reads per sample in the model. Here, same thing as the PERMANOVA: replace Treatment and Time with the factors you want to test, leaving what you are more interested in as first factor (e.g. if you are more interested in the Treatment rather than the Time, your model will be Treatment + Time. If you want to know also their interaction your model will be Treatment * Time). The order of your factors is particular important when you have an unbalanced design.

rowSums(ASV_data_cleaned)
fact$numberReads <- rowSums(ASV_data_cleaned)

#use log number of reads for quicker result
dat_nb_spat <- manyglm(dat_mvabund ~ Treatment * Family + offset(log(numberReads)), family = "negative.binomial", data = fact)

dat_ps <- manyglm(dat_mvabund ~ Treatment * Age * Family + offset(log(numberReads)), family = "poisson", data = fact)

# 
# Option with COMPOSITION = TRUE 
# This is the best option when dealing with compositional data!!!! It runs slower but we should always use COMPOSITION = TRUE. 
# This replaces the log(numberReads)

dat_nb_compositionT = manyglm(dat_mvabund ~ Treatment * Age * Family, family="negative.binomial", data = fact, composition=TRUE)


#If running into Error: vector memory exhausted (limit reached?)
#To solve go into your computer terminal and put "cd ~ touch .Renviron open .Renviron"
#Then in new window put "R_MAX_VSIZE=100Gb"
#This changes computer settings to allocate more space to Rstudio

# # check model assumptions
# # check the mean-variance relationship
pdf("mb2021_plots/meanvar_plot.pdf", width = 5, height = 5)
meanvar.plot(dat_nb, legend = TRUE)
# add line for Poisson mean-variance relationship
abline(0,1)
dev.off()
# 

# # Check the model fitting plotting the residuals
plot(dat_nb, n.vars = 50) # involves some randomisation, repeat it a couple of times
# if there are No obvious pattern among multiple plots - it is a good model fit.
plot(dat_nb, which = 2)
# 
# # univariate models # This is a VERY LONG STEP!!! I always run it in Katana. 
# You can decide to use "unadjusted" or "adjusted" p-values. I usually run it with both ans see the differences later. Somethimes adjusted can be too strict and give me only very few significant zASV, on the other hand, unadjusted gives me too many. So I run both an see what is the best for my  data.

pseq <- mb2021_filtered_NOT_rarefied

pseq <- subset_samples(pseq, Age %in% c("1 dpf"))

#correct family names
pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

fact1 = sample_data(pseq)
fact = as.matrix.data.frame(fact1)
fact = as.data.frame(fact)
ASV_data_cleaned <- pseq@otu_table

ASV_data_cleaned <- as.data.frame(pseq@otu_table)
rowSums(ASV_data_cleaned)
fact$numberReads <- rowSums(ASV_data_cleaned)
dat_mvabund <- mvabund(ASV_data_cleaned)

dat_nb_day1 <- manyglm(dat_mvabund ~ Treatment * Family + offset(log(numberReads)), family = "negative.binomial", data = fact)

dat_aov_unadj_day1 <- anova(dat_nb_day1, p.uni = "unadjusted") # This is a VERY LONG STEP!!! I always run it in Katana. 

#option 2
dat_aov_adj_spat <- anova(dat_nb_spat, p.uni = "adjusted") # This is a VERY LONG STEP!!! I always run it in Katana. 
dat_aov_adj$table

length(which(dat_aov_adj_spat$uni.p[2,] < 0.05)) # Number of ASVs affected by treatment
length(which(dat_aov_adj_spat$uni.p[3,] < 0.05)) # Number of ASVs affected by time
length(which(dat_aov_adj_spat$uni.p[4,] < 0.05)) # Number of ASVs affected by interaction treatment * time

save.image("data.RData") #This is IMPORTANT! If you save this you won't have to run the long step again!!!

summary(dat_nb) #long step too.

diff_05_treat_unj_spat <- which(dat_aov_unadj_spat$uni.p[2,] < 0.05) # [2] is referring to the first factor of my model , in my case "Treatment"
diff_05_time_unj <- which(dat_aov_unadj$uni.p[3,] < 0.05) # [3] is referring to the second factor of my model , in my case "Time"
diff_05_interact_unj <- which(dat_aov_unadj$uni.p[4,] < 0.05) # [4] is referring to the third factor of my model , in my case "interaction Treatment and Time"
# 
# #### Compare ASVs (in your case will be the number of molecular features) affected by treatment, time and those affected by the interaction ####
names_diff_treat_unj_spat <- names(diff_05_treat_unj_spat)
names_diff_time_unj <- names(diff_05_time_unj)
names_diff_int_unj <- names(diff_05_interact_unj)
# 
length(which(names_diff_treat %in% names_diff_time))

***stopped here Friday****

# # Number of ASVs affected by both treatment and time
round(length(which(names_diff_treat %in% names_diff_time))/ncol(ASV_data_cleaned) * 100) # % of ASVs affected by both treatment and time
# 
length(which(names_diff_treat %in% names_diff_int))
# #  Number of ASVs affected by both treatment and interaction
length(which(names_diff_time %in% names_diff_int))
# # Number of ASVs affected by both time and interaction
# 
# #### saving outputs -- ASV counts -- raw data ####
dat_diff_treat <- ASV_data_cleaned[, which(names(ASV_data_cleaned) %in% names_diff_treat)]
dim(dat_diff_treat) # Number of ASVs affected by treatment
dat_diff_treat$response <- "treatment"
write.csv(dat_diff_treat, "outputs/Files/rawData_diffabundantASVs_mvabund_treatment.csv")
round(dim(dat_diff_treat)[2]/ncol(ASV_data_cleaned) * 100) # % ASVs affected by treatment
# 
dat_diff_time <- ASV_data_cleaned[, which(names(ASV_data_cleaned) %in% names_diff_time)]
dim(dat_diff_time) # Number ASVs affected by Time
dat_diff_time$response <- "time"
# write.csv(dat_diff_time, "outputs/Files/rawData_diffabundantASVs_mvabund_time.csv")
round(dim(dat_diff_time)[2]/ncol(ASV_data_cleaned) * 100) # % ASVs affected by Time
# 
dat_diff_int <- ASV_data_cleaned[, which(names(ASV_data_cleaned) %in% names_diff_int)]
dim(dat_diff_int) # Number ASVs affected by interaction
dat_diff_int$response <- "interactionTreatmentTime"
# write.csv(dat_diff_int, "outputs/Files/rawData_diffabundantASVs_mvabund_interactionTreatmentTime.csv")
round(dim(dat_diff_int)[2]/ncol(ASV_data_cleaned) * 100) # % ASVs affected by interaction
# 
# # find ASVs that are responsive to more than one condition (treatment, time or interaction)
# # combine datasets
dat_diff_treat_t <- as.data.frame(dat_diff_treat %>%
                                    dplyr::select(-response) %>%
                                    t()) %>%
  mutate(response = "treatment", ASV_ID = rownames(.))

dat_diff_time_t <- as.data.frame(dat_diff_time %>%
                                   dplyr::select(-response) %>%
                                   t()) %>%
  mutate(response = "time", ASV_ID = rownames(.))
dat_diff_int_t <- as.data.frame(dat_diff_int %>%
                                  dplyr::select(-response) %>%
                                  t()) %>%
  mutate(response = "interactionTreatmentTime", ASV_ID = rownames(.))

all_diff_ASVs <- rbind(dat_diff_treat_t, dat_diff_time_t, dat_diff_int_t)
head(all_diff_ASVs)
tail(all_diff_ASVs)

all_diff_ASVs <- all_diff_ASVs %>%
  dplyr::select(ASV_ID, everything())

length(unique(all_diff_ASVs$ASV_ID)) # Number of unique ASVs
round(length(unique(all_diff_ASVs$ASV_ID))/ncol(ASV_data_cleaned) * 100) # % of unique ASVs
# 
# # Find unique ASVs - those affected by single factor - either treatment, time OR interaction
unique_diff_ASVs <- all_diff_ASVs %>%
  group_by(ASV_ID) %>%
  filter(n() == 1) %>%
  as.data.frame()
# 
dim(unique_diff_ASVs) # 3267   56
length(unique(unique_diff_ASVs$ASV_ID)) # 3267

dim(unique_diff_ASVs %>% filter(response == "treatment")) # Number ASV affected only by treatment
dim(unique_diff_ASVs %>% filter(response == "time")) #  Number ASV affected only by time
dim(unique_diff_ASVs %>% filter(response == "interactionTreatmentTime")) #  Number ASV affected only by interaction

unique_diff_ASVs$combined_response <- unique_diff_ASVs$response
head(unique_diff_ASVs)
# 
# # find duplicated ASVs - those affected in more than one scenario (treatment, time, interaction)
duplic_diff_ASVs <- as.data.frame(all_diff_ASVs %>%
                                    group_by(ASV_ID) %>%
                                    filter(n() > 1))
dim(duplic_diff_ASVs) # 
length(unique(duplic_diff_ASVs$ASV_ID)) # Number unique ASVs
# 
# # find ASVs affected in the 3 scenarios (treatment, time AND interaction)
duplic_diff_ASVs_ALL <- as.data.frame(duplic_diff_ASVs %>%
                                        group_by(ASV_ID) %>%
                                        filter(n() == 3)) %>%
  arrange(ASV_ID)

head(duplic_diff_ASVs_ALL)
dim(duplic_diff_ASVs_ALL) # 
length(unique(duplic_diff_ASVs_ALL$ASV_ID)) # Number ASVs affected by 3 scenarios
# 
duplic_diff_ASVs_ALL$combined_response <- "treatment_time_interactionTreatmentTime"
# 
# # find ASVs affected in by treatment and time
duplic_ASVs_treatTime <- as.data.frame(duplic_diff_ASVs %>%
                                         group_by(ASV_ID) %>%
                                         filter(n() == 2, response == "treatment" | response == "time") %>%
                                         ungroup() %>%
                                         group_by(ASV_ID) %>%
                                         filter(n() == 2) %>%
                                         ungroup() %>%
                                         arrange(ASV_ID))
# 
head(duplic_ASVs_treatTime)
dim(duplic_ASVs_treatTime) # 
length(unique(duplic_ASVs_treatTime$ASV_ID)) # Number of ASVs affected by treatment and time

duplic_ASVs_treatTime$combined_response <- "treatment_time"
round(length(unique(duplic_ASVs_treatTime$ASV_ID))/ ncol(ASV_data_cleaned) * 100) # %
# 
# # find ASVs affected in by treatment and interaction
duplic_ASVs_treatInter <- as.data.frame(duplic_diff_ASVs %>%
                                          group_by(ASV_ID) %>%
                                          filter(n() == 2, response == "treatment" | response == "interactionTreatmentTime") %>%
                                          ungroup() %>%
                                          group_by(ASV_ID) %>%
                                          filter(n() == 2) %>%
                                          ungroup() %>%
                                          arrange(ASV_ID))
# 
head(duplic_ASVs_treatInter)
dim(duplic_ASVs_treatInter) # 
length(unique(duplic_ASVs_treatInter$ASV_ID)) # Number ASVs affected by treatment and interactionTreatmentTime

duplic_ASVs_treatInter$combined_response <- "treatment_interactionTreatmentTime"

# # find ASVs affected in by time and interaction
duplic_ASVs_timeInter <- as.data.frame(duplic_diff_ASVs %>%
                                         group_by(ASV_ID) %>%
                                         filter(n() == 2, response == "time" | response == "interactionTreatmentTime") %>%
                                         ungroup() %>%
                                         group_by(ASV_ID) %>%
                                         filter(n() == 2) %>%
                                         ungroup() %>%
                                         arrange(ASV_ID))

head(duplic_ASVs_timeInter)
dim(duplic_ASVs_timeInter) # 
length(unique(duplic_ASVs_timeInter$ASV_ID)) # Number ASVs affected by time and interactionTreatmentTime

duplic_ASVs_timeInter$combined_response <- "time_interactionTreatmentTime"

#make sure all objects for rbind are populatated (aka have signif ASVs and are not equal to zero)
diff_ASVs_final <- rbind(unique_diff_ASVs, duplic_diff_ASVs_ALL, duplic_ASVs_treatTime, duplic_ASVs_treatInter, duplic_ASVs_timeInter)

dim(diff_ASVs_final) # 
diff_ASVs_final <- diff_ASVs_final %>%
  arrange(ASV_ID) %>%
  dplyr::select(-response)
head(diff_ASVs_final)

write.csv(diff_ASVs_final, "outputs/files/diffabundantASVs_int_unadj_mvabund.csv", row.names = FALSE)
# 


#Normalize data for figures ---- 
#Data normalisation using DESeq2:
#counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") #download deseq

BiocManager::install("DESeq2") #Install deseq

BiocManager::install("DESeq2", force = TRUE)

BiocManager::install("GenomeInfoDb")
BiocManager::install("DESeq2")

library(DESeq2)
library(phyloseq)

mb2021_filtered_NOT_rarefied <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/mb2021_filtered_NOT_rarefied.rds")

pseq <- mb2021_filtered_NOT_rarefied

DeSeq <- phyloseq_to_deseq2(pseq, ~ Treatment) #convert phyloseq to deseq object

DeSeq2 <- DESeq(DeSeq)

# # two-way analysis with 2 levels each (treatment + time)
# # countData = ASV counts with samples in columns  
# # colData = Sample meta data (factors, etc)  
# # design = the right hand side of a GLM formula

ASV_deseq <- DESeqDataSetFromMatrix(countData = t(),
                                    colData = fact2,
                                    design = ~ Treatment)

# positive counts normalisaton (accounts for zero-inflation)
ASV_deseq <- estimateSizeFactors(DeSeq2, type = "poscounts")
sizeFactors(ASV_deseq)

# size-factor corrected data are calculated by dividing the raw counts by the sample size factor and adding 0.5 to correct for the zeros
ASV_deseq_norm <- sapply(row.names(ASV_deseq), function(x){
  plotCounts(ASV_deseq, x, "Treatment", returnData = TRUE, normalized = TRUE)$count
  
})
rownames(ASV_deseq_norm) <- colnames(ASV_deseq)

# remove the addition of 0.5 to all entries
ASV_deseq_norm <- ASV_deseq_norm - 0.5
ASV_deseq_norm[1:10, 1:10]

# make dataset with integers
ASV_deseq_norm <- round(ASV_deseq_norm)
ASV_deseq_norm[1:10, 1:10]

View(ASV_deseq_norm)

# check for zeros across all samples (produced by normalisation)
dim(ASV_deseq_norm[, colSums(ASV_deseq_norm) == 0]) # 0

# check singletons
dim(ASV_deseq_norm[, colSums(ASV_deseq_norm) == 1])

#save normalised table
write.csv(ASV_deseq_norm, "~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files//normalised_ASV_table.csv")

#Make new pseq object with normalized ASV table
mb2021_filtered_NOT_rarefied <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/mb2021_filtered_NOT_rarefied.rds")

#Check format
View(mb2021_filtered_NOT_rarefied@otu_table)
View(normalised_ASV_table)

#correct to make forst column as row names

normalised_ASV_table <- as.data.frame(normalised_ASV_table)

rownames(normalised_ASV_table) <- normalised_ASV_table[, 1]

# Remove the first column from the data frame
normalised_ASV_table <- normalised_ASV_table[, -1]

# Create a new otu_table object
new_otu_table <- otu_table(normalised_ASV_table, taxa_are_rows = FALSE)

# Replace the old OTU table with the new one
otu_table(mb2021_filtered_NOT_rarefied) <- new_otu_table

# Verify the replacement
View(mb2021_filtered_NOT_rarefied@otu_table)
#Save as different RDS file
saveRDS(mb2021_filtered_NOT_rarefied, file = "mb2021_filtered_NOT_rarefied_normalized.rds")

