#ANCOM-BC 
#https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html 

#ANCOM-BC is used as a bias correction to control for sampling fraction and sequencing efficiency biases
#Ancom-BC performs a log transformation and runs a sensitivity analyses for the value used for a pseudo-count (on zeros)

#Performing ancom on data that has been cleaned (e.g., low prev taxa remove, cholor/mito removed) but
#sequence depth has not been rarefied AND not deseq2 normalized 

#Packages

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ANCOMBC")
# install.packages("DT")

library(ANCOMBC)
packageVersion("ANCOMBC")

library(tidyverse)
library(DT)


#Something is wrong with how I am using ANCOMBC....
#NOTHING comes up as signficantly different, even when I test for differences
#between age groups. 
#And sometime the lfc is completely wrong and does not make sense. 
#Wondering if my data does not work well as high zero counts in the data

pseq <- PB2023_spat_ANCOMBC #Nothing significant
pseq <- PB2023_spat_not_rarefied_normalized #Nothing significant
pseq <- PB2023_spat_limited_10X_reads
pseq <- MU42022_normalized_relative #deseq normalized

pseq <- MU42022_filtered_Oct92024
pseq <- PB2023_spat_filtered_not_rarefied
pseq <- PB2023_rarefied_3000

pseq <- mb2021_filtered_NOT_rarefied

# pseq@sam_data$Treatment <- gsub("Killed-Probiotics", "KilledProbiotics", pseq@sam_data$Treatment)
# pseq@sam_data$Treatment

# Subset the dataset to include only the relevant treatments
pseq <- subset_samples(pseq, Treatment %in% c("Control", "Probiotics", "Killed-Probiotics"))

#MU42022 filtering
# pseq <- subset_samples(pseq, Age %in% c("Spat"))
# pseq <- subset_samples(pseq, !Genetics %in% c("4"))
# pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))
# 

#Sanity check
#check if any OTUs are not present in any samples (want false)
any(taxa_sums(pseq) == 0)

#if true

pseq_filtered <- prune_taxa(taxa_sums(pseq) > 0, pseq)
any(taxa_sums(pseq_filtered) == 0)

pseq <- pseq_filtered

#Permutation
set.seed(123)

pseq_perm = pseq
meta_data_perm = microbiome::meta(pseq_perm)
meta_data_perm$Treatment = sample(meta_data_perm$Treatment)
meta_data_perm$Family <- sample(meta_data_perm$Family)
phyloseq::sample_data(pseq_perm) = meta_data_perm

#All factors
# meta_data_perm <- microbiome::meta(pseq_perm)
# meta_data_perm$Treatment <- sample(meta_data_perm$Treatment)
# meta_data_perm$Age <- sample(meta_data_perm$Age)
# # meta_data_perm$Genetics <- sample(meta_data_perm$Genetics)
# meta_data_perm$Family <- sample(meta_data_perm$Family)
# phyloseq::sample_data(pseq_perm) <- meta_data_perm


#Run ancombc fxn

output = ancombc2(data = pseq_perm, tax_level = "Genus",
                  fix_formula = "Treatment",
                  p_adj_method = "holm", 
                  pseudo_sens = TRUE,
                  prv_cut = 0.2, #ASVs present in < 20% removed
                  lib_cut = 1000, 
                  s0_perc = 0.05,
                  group = "Treatment",
                  struc_zero = TRUE, 
                  neg_lb = TRUE,
                  dunnet = TRUE)

#Random effect - not converging
output = ancombc2(data = pseq_perm, tax_level = "Genus",
                  fix_formula = "Treatment", 
                  rand_formula = "Treatment | Family",
                  p_adj_method = "bonferroni", 
                  pseudo_sens = TRUE,
                  prv_cut = 0.1, 
                  lib_cut = 1000, 
                  s0_perc = 0.05,
                 group = "Treatment",
                 lme_control = lme4::lmerControl(),
                 struc_zero = FALSE, 
                 neg_lb = FALSE)


# Run ANCOMBC2 with all factors
output <- ancombc2(
  data = pseq_perm,
  tax_level = "Genus",
  fix_formula = "Treatment + Age + Genetics", # Include all desired variables
  p_adj_method = "holm",
  pseudo_sens = FALSE,
  prv_cut = 0,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = NULL, # Set to NULL if not grouping by one variable
  struc_zero = FALSE,
  neg_lb = FALSE
)


#Structural zeros table (presence/absence)
tab_zero = output$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")


#primary output of the ANCOM-BC2 methodology identifies taxa with differential 
#abundance based on the chosen covariate. 
#The results include: 1) log fold changes, 2) standard errors, 3) test statistics, 
#4) p-values, 5) adjusted p-values, 6) indicators denoting whether the taxon is 
#differentially abundant (TRUE) or not (FALSE), and 7) indicators denoting whether the 
#taxon passed the sensitivity analysis (TRUE) or not (FALSE)

res_prim = output$res
View(res_prim)

# Filter for significant genera
significant_genera <- res_prim %>%
  filter(p_TreatmentProbiotics < 0.05 | 'p_TreatmentProbiotics+HT' < 0.05)
print(significant_genera)

library(ggplot2)

# Filter for Celeribacter
celeribacter_data <- res_prim[res_prim$taxon == "Phaeobacter", ]

# Create a tidy data frame for plotting
lfc_data <- data.frame(
  Treatment = c("Control", "Probiotics", "Probiotics + HT"),
  LFC = c(
    celeribacter_data$`lfc_(Intercept)`,  # Note the backticks for special characters
    celeribacter_data$`lfc_TreatmentProbiotics`,
    celeribacter_data$`lfc_TreatmentProbiotics + HT`
  )
)

# Check the data
print(lfc_data)

# Plot the LFC values
ggplot(lfc_data, aes(x = Treatment, y = LFC, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  labs(
    title = "Log-Fold Change (LFC) of Celeribacter",
    x = "Treatment",
    y = "Log-Fold Change (LFC)"
  )


#Something is wrong - says there is lower PB (ASV7/18) in PB treatments





