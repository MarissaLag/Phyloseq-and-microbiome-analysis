#Zicoseq differential abundace test
#see paper: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01320-0
#Supposed to control for FDR but keep higher power (unlike ANCOMBC and ALDEx2)
#See package info: https://cran.r-project.org/web/packages/GUniFrac/vignettes/ZicoSeq.html

#load
install.packages("GUniFrac")
library(GUniFrac)

#Set theme
theme.marissa <- function() {
  theme_classic(base_size = 14) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"))
}

theme_set(theme.marissa())

#pseq object - data

#may want to use CSS normalized 
pseq <- PB2023_spat_not_rarefied_CSSnormalized_Jan2025
pseq <- PB2023_spat_filtered_not_rarefied
pseq <- PB2023_rarefied_3000
pseq <- PB2023_spat_10X_limited_CSS
pseq <- MU42022_filtered_Oct92024
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Day %in% c("Spat"))
pseq <- subset_samples(pseq, !Organism %in% c("Water"))
pseq <- subset_samples(pseq, Microbial.Source %in% c("Disease Resilient", "Disease Susceptible", "Control"))
pseq <- subset_samples(pseq, !Treatment %in% c("James", "Continuous Probiotics"))

pseq <- MU42022_filtered_NOT_rarefied
pseq <- mb2021_filtered_NOT_rarefied_normalized
pseq <- mb2021_filtered_NOT_rarefied_June2024
pseq <- mb2021_filteredw1dpf_only_rarefied_June2024
pseq <- Sam_spatonly_Rare
pseq <- subset_samples(pseq, Age %in% c("Spat"))

pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

pseq <- subset_samples(pseq, !Family %in% c("1"))


#Quality check
any(taxa_sums(pseq) == 0)

#if true

pseq_filtered <- prune_taxa(taxa_sums(pseq) > 0, pseq)
any(taxa_sums(pseq_filtered) == 0)

pseq <- pseq_filtered

OTU <- pseq@otu_table
Meta <- pseq@sam_data

comm <- t(OTU)

meta.dat = sample_data(pseq)
meta = as.matrix.data.frame(meta.dat)
meta.dat = as.data.frame(meta)
str(meta.dat)

meta.dat$Treatment <- as.factor(meta.dat$Treatment) 
meta.dat$Family <- as.factor(meta.dat$Family)
meta.dat$Genetics <- as.factor(meta.dat$Genetics)

meta.dat$Genetic.Background <- as.factor(meta.dat$Genetic.Background) 
meta.dat$Microbial.Source <- as.factor(meta.dat$Microbial.Source)

meta.dat$Microbial.Source <- relevel(factor(meta.dat$Microbial.Source), ref = "Control")

#RUN TEST
ZicoSeq.obj <- ZicoSeq(meta.dat = meta.dat, feature.dat = comm, 
                       grp.name = 'Microbial.Source', feature.dat.type = "other",
                       adj.name = "Genetic.Background",
                       # Filter to remove rare taxa
                       prev.filter = 0.1, mean.abund.filter = 0,  
                       max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.01, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 3, 
                       # Use the square-root transformation - recommended by author
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 9999,  strata = NULL, 
                       # Reference-based multiple stage normalization
                       ref.pct = 0.3, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)


#For sample size less than 40, posterior sampling will not be used!
# 53  features are filtered!
#   The data has  17  samples and  130  features will be tested!
#   On average,  1  outlier counts will be replaced for each feature!

names(ZicoSeq.obj)

# Raw p-values
raw_p_values <- ZicoSeq.obj$p.raw
significant_asvs_raw <- names(raw_p_values)[raw_p_values < 0.01]

# Adjusted p-values for FDR
adj_p_values_fdr <- ZicoSeq.obj$p.adj.fdr
significant_asvs_fdr <- names(adj_p_values_fdr)[adj_p_values_fdr < 0.05]

# Adjusted p-values for FWER
adj_p_values_fwer <- ZicoSeq.obj$p.adj.fwer
significant_asvs_fwer <- names(adj_p_values_fwer)[adj_p_values_fwer < 0.05]


tapply(rowSums(comm), meta.dat$Treatment, mean)

#plots
ZicoSeq.plot(ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = 0.05, text.size = 10,
             out.dir = NULL, width = 10, height = 6) 

ZicoSeq.plot(ZicoSeq.obj, pvalue.type = 'p.raw', cutoff = 0.05, text.size = 10,
             out.dir = NULL, width = 10, height = 6)

ZicoSeq.plot(ZicoSeq.obj, pvalue.type = 'p.adj.fwer', cutoff = 0.05, text.size = 10,
             out.dir = NULL, width = 10, height = 6)


#Editing zicoSeq to see if FDR can be less strict
ZicoSeq.obj <- ZicoSeq(meta.dat = meta.dat, feature.dat = comm, 
                       grp.name = 'Treatment', feature.dat.type = "count",
                       # Filter to remove rare taxa
                       prev.filter = 0.2, max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = FALSE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 9999, strata = NULL, 
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = FALSE, verbose = TRUE, return.feature.dat = TRUE)

print(length(ZicoSeq.obj$p.raw))  # Number of ASVs tested - 130, so ~50 removed

#Disable normalization
ZicoSeq.obj <- ZicoSeq(meta.dat = meta.dat, feature.dat = comm, 
                       grp.name = 'Treatment', feature.dat.type = "count",
                       prev.filter = 0.1, max.abund.filter = 0.001, min.prop = 0,  
                       is.winsor = TRUE, outlier.pct = 0.05, winsor.end = 'top',
                       is.post.sample = FALSE, post.sample.no = 25, 
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       perm.no = 9999, strata = NULL, 
                       ref.pct = 0, stage.no = 1,  # 🔹 Disable normalization
                       is.fwer = FALSE, verbose = TRUE, return.feature.dat = TRUE)

#Look at abundance of signif ASVs

# Load necessary libraries
library(phyloseq)
library(ggplot2)

pseq <- mb2021_filtered_NOT_rarefied_normalized
pseq <- subset_samples(pseq, Age %in% c("Spat"))

pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

# Define the ASVs of interest

significant_asvs <- significant_asvs_fdr

significant_asvs <- c("ASV190",  "ASV227",  "ASV231",  "ASV236",
                      "ASV348",  "ASV461",  "ASV478",  "ASV777",  "ASV953",  "ASV1211")

significant_asvs <- c("ASV6", "ASV20", "ASV41", "ASV42", "ASV77",
                      "ASV101", "ASV126", "ASV136", "ASV249", "ASV263", 
                      "ASV267", "ASV276", "ASV282", "ASV322", "ASV346", "ASV413")

significant_asvs <- c("ASV60", "ASV101")

# Subset the phyloseq object to include only the ASVs of interest
pseq_subset <- prune_taxa(taxa_names(pseq) %in% significant_asvs, pseq)

# If necessary, transform the data (e.g., relative abundance)
pseq_subset <- transform_sample_counts(pseq_subset, function(x) x / sum(x))

# Convert to data frame for ggplot
df <- psmelt(pseq_subset)

# Plotting the abundance of significant ASVs across treatments
ggplot(df, aes(x = Microbial.Source, y = Abundance, fill = Microbial.Source)) + # or another relevant taxonomy level
  geom_boxplot() +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) +  # Custom color for Treatment
  facet_wrap(~ OTU, scales = "free") + # creates a facet for each ASV
  labs(title = "Abundance of Celeribacter - Spat",
       x = "",
       y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

