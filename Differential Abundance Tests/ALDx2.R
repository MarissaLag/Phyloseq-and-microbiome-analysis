#USing ALDEx2 for differntial abundance analysis
#Indic species uses count-based rather than compositional methods
#ALDEx2 uses compositional methods

#ALDEx2 expects count data (not compositional) - so don't transform into compositional first

install.packages("ALDEx2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ALDEx2")

library(ALDEx2)
library(phyloseq)

pseq <- MU42022_filtered_Oct92024

pseq <- PB2023_spat_not_rarefied_normalized

pseq <- MU42022_nonrare_normalized_spat_limitedtaxa

pseq <- PB2023_spat_not_rarefied_CSSnormalized_Jan2025

pseq <- PB2023_spat_filtered_not_rarefied

pseq <- PB2023_rarefied_3000

pseq <- PB2023_spat_limited_10X_reads

#MU42022 filtering
# pseq <- subset_samples(pseq, Age %in% c("Spat"))
# pseq <- subset_samples(pseq, !Genetics %in% c("4"))
# pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))

#PB2023 filtering
pseq <- subset_samples(pseq, !Treatment %in% c("James", "Continuous Probiotics"))


# Extract the OTU table and convert it to a matrix
otu_table_matrix <- as(otu_table(pseq), "matrix")

#If CSS normalized, must covert to integers
# Rounding the values to the nearest integer
otu_table_matrix <- round(otu_table_matrix)

# Transpose the OTU table matrix
otu_table_matrix <- t(otu_table_matrix)

# Extract sample data to get the group information
sample_info <- sample_data(pseq)

sample_info$Treatment <- as.factor(sample_info$Treatment)
group_vector <- sample_info$Treatment 
group_vector <- as.character(group_vector)


# Perform CLR transformation and differential abundance test with ALDEx2
# 'mc.samples' is the number of Monte Carlo instances, typically set to 128 or 500
#Use 500 for smaller sample sizes
aldex_clr <- aldex.clr(otu_table_matrix, group_vector, mc.samples = 500)

# Use ALDEx2 for a t-test or ANOVA-like test depending on your groups
#t-test is for comparing 2 factors, Anova for more than 2 factors:
# For two groups: use aldex.ttest; for multiple groups: use aldex.kw
aldex_results <- aldex.kw(aldex_clr)
head(aldex_results)

# Filter for significant ASVs based on eBH values < 0.05
significant_kw <- aldex_results[aldex_results$kw.ep < 0.05, ]
significant_glm <- aldex_results[aldex_results$glm.ep < 0.05, ]

# Optionally, you can combine these results 
significant_all <- merge(significant_kw, significant_glm, by = "row.names", all = TRUE)

# View results
head(significant_kw)
head(significant_glm)
head(significant_all)

# Example boxplot for one ASV
significant_asv <- c("ASV190", "ASV227")  # Replace with actual significant ASV name

# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(dplyr)

# Assuming 'ps' is your phyloseq object
# Melt the phyloseq object to a data frame
pseq <- subset_samples(pseq, !Treatment %in% c("Continuous Probiotics", "James"))
ps_melted <- psmelt(pseq)

# Filter for ASV9
asv9_data <- ps_melted %>% filter(OTU == significant_asv)

# Calculate mean and standard error for each treatment
summary_data <- asv9_data %>%
  group_by(Treatment, OTU) %>%
  summarise(
    mean_abundance = mean(Abundance),
    se_abundance = sd(Abundance) / sqrt(n())  # Standard Error
  )

# Create the plot
ggplot(summary_data, aes(x = Treatment, y = mean_abundance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
                width = 0.2) +
  labs(title = "",
       x = "Treatment",
       y = "Mean Abundance (± SE)") +
  facet_wrap(~OTU) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#Zero-inflation may be effecting test?
#Use 
otu_table <- pseq@otu_table
# Filter out ASVs that are present in fewer than 10% of samples
otu_table_filtered <- otu_table[rowSums(otu_table > 0) / ncol(otu_table) > 0.1, ]


# Assuming 'otu_table' is your count data (samples in rows, ASVs in columns)
otu_table <- as.data.frame(pseq@otu_table)
View(otu_table)
str(otu_table)

library(ggplot2)
library(reshape2)

# Create a heatmap of zero counts
zero_matrix <- otu_table == 0
zero_melted <- melt(zero_matrix)
ggplot(zero_melted, aes(Var1, Var2, fill=value)) +
  geom_tile() + 
  theme_minimal() + 
  scale_fill_manual(values=c("white", "black"))

# Choose a small value to add to zeros
small_value <- 1

# Replace all zeros in the otu_table with the small_value
otu_table[otu_table == 0] <- small_value

# Check the updated data
head(otu_table)

otu_table_matrix <- t(otu_table)

#try ALD again
# Perform CLR transformation and differential abundance test with ALDEx2
# 'mc.samples' is the number of Monte Carlo instances, typically set to 128 or 500
aldex_clr <- aldex.clr(otu_table_matrix, group_vector, mc.samples = 128, )

# Use ALDEx2 for a t-test or ANOVA-like test depending on your groups
#t-test is for comparing 2 factors, Anova for more than 2 factors:
# For two groups: use aldex.ttest; for multiple groups: use aldex.kw
# Here we assume two groups (use aldex.kw if more than two)
aldex_results <- aldex.kw(aldex_clr)

head(aldex_results)

# Filter for significant ASVs based on eBH values < 0.05
significant_kw <- aldex_results[aldex_results$kw.eBH < 0.05, ]
significant_glm <- aldex_results[aldex_results$glm.eBH < 0.05, ]

# View results
head(significant_kw)
head(significant_glm)
head(significant_all)

significant_asvs <- rownames(aldex_results[aldex_results$glm.ep < 0.05, ])  # ASVs with adjusted p-value < 0.05

#Pairwise analysis
#Can only handle 2 groups
x.effect <- aldex.effect(aldex_clr, CI=T, verbose=F, include.sample.summary=T, 
                         paired.test=T)
