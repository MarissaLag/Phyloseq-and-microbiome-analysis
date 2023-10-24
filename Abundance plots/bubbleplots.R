#Look at pseq.core

#Quick lib ----

library("devtools")
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)


#Load data ----

Marissa_MU42022_rarefied_20231016 <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rarefied_20231016.rds")

pseq <-  Marissa_MU42022_rarefied_20231016

#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#Day 1 only ----

pseq <- subset_samples(pseq, !Age %in% c("Spat", "Day 03", "Day 06", "Day 15"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#Spat only ----

pseq <- subset_samples(pseq, !Age %in% c("Day 01", "Day 03", "Day 06", "Day 15"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#remove algae ----

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#Larvae only with algae ----

pseq <- subset_samples(pseq, !Age %in% "Spat")

#Larvae only without algae ----

pseq <- subset_samples(pseq, !Age %in% "Spat")

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")


#Top phyla, all samples ----

top5P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Phylum"], sum), TRUE)[1:5]

top5P = subset_taxa(pseq, Phylum %in% names(top5P.names))

pseq <- top5P

#Create pseq objects ----

pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq, level = "Genus", detection = 50/100, prevalence = 70/100)

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 95/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")


#pseq core - Probiotics + PB + Heat only (must remove all other samples first

pseq <- subset_samples(pseq, !Treatment %in% c("Control", "High temperature"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

pseq.core <- core(pseq, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

#plot MDS/PcoA ----

set.seed(4235421)

ord <- ordinate(pseq, "MDS", "bray")

plot_ordination(pseq, ord, color = "Treatment", shape = "Age") + geom_point(size = 4)


#psmelt ----

#all data

pseq <- psmelt(pseq)

pseq_core <- psmelt(pseq.core)

pseq_fam <- psmelt(pseq_fam)

pseq_phy <- psmelt(pseq_phy)

pseq_phy_rel <- psmelt(pseq.phy.rel)

pseq_fam_rel <- psmelt(pseq.fam.rel)

pseq_gen_rel <- psmelt(pseq.gen.rel)


#Day 1

pseq_Day1 <- psmelt(pseq)

pseq_core_Day1 <- psmelt(pseq.core)

pseq_fam_Day1 <- psmelt(pseq_fam)

pseq_phy_Day1 <- psmelt(pseq_phy)

pseq_phy_rel_Day1 <- psmelt(pseq.phy.rel)

pseq_fam_rel_Day1 <- psmelt(pseq.fam.rel)

pseq_gen_rel_Day1 <- psmelt(pseq.gen.rel)

#spat

pseq_spat <- psmelt(pseq)

pseq_core_spat <- psmelt(pseq.core)

pseq_fam_spat <- psmelt(pseq_fam)

pseq_phy_spat <- psmelt(pseq_phy)

pseq_phy_rel_spat <- psmelt(pseq.phy.rel)

pseq_fam_rel_spat <- psmelt(pseq.fam.rel)

pseq_gen_rel_spat <- psmelt(pseq.gen.rel)

#remove algae

pseq_oyster <- psmelt(pseq)

pseq_core_oyster <- psmelt(pseq.core)

pseq_fam_oyster <- psmelt(pseq_fam)

pseq_phy_oyster <- psmelt(pseq_phy)

pseq_phy_rel_oyster <- psmelt(pseq.phy.rel)

pseq_fam_rel_oyster <- psmelt(pseq.fam.rel)

pseq_gen_rel_oyster <- psmelt(pseq.gen.rel)

#larvae only with algae

pseq_larvae_a <- psmelt(pseq)

pseq_core_larvae_a <- psmelt(pseq.core)

pseq_fam_larvae_a <- psmelt(pseq_fam)

pseq_phy_larvae_a <- psmelt(pseq_phy)

pseq_phy_rel_larvae_a <- psmelt(pseq.phy.rel)

pseq_fam_rel_larvae_a <- psmelt(pseq.fam.rel)

pseq_gen_rel_larvae_a <- psmelt(pseq.gen.rel)

#larvae only without algae

pseq_larvae <- psmelt(pseq)

pseq_core_larvae <- psmelt(pseq.core)

pseq_fam_larvae <- psmelt(pseq_fam)

pseq_phy_larvae <- psmelt(pseq_phy)

pseq_phy_rel_larvae <- psmelt(pseq.phy.rel)

pseq_fam_rel_larvae <- psmelt(pseq.fam.rel)

pseq_gen_rel_larvae <- psmelt(pseq.gen.rel)


#Top phyla

Top5P <- psmelt(top5P)

#If not using psmelt: 
#Extract core data as .csv ----


file_path <- "C:/Users/maris/OneDrive/Documents/GitHub/pseq_core_spat.csv"

# Write the otu_table to a CSV file
write.csv(pseq.core@otu_table, file = file_path)

if (file.exists(file_path)) {
  cat("OTU table CSV file has been created:", file_path, "\n")
} else {
  cat("Failed to create the CSV file.")
}

#Have to manually change sample names into treatments in excel sheet
#import excel sheet - have to make it with x, y, z setup (data as rows)





###Bubbble plots ----


###pseq ----

ggplot(pseq_core, aes(x=Treatment, y=OTU, size = Abundance, color = Family)) + 
  geom_point(alpha=0.7)+ 
  scale_size(range = c(.3, 12)) +
  scale_colour_ipsum() +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("Taxa (Family)") +
  xlab("") +
  theme(legend.position = "none")

###pseq relative ----

ggplot(pseq_fam_rel_spat, aes(x=Treatment, y=Family, size = Abundance, color = Treatment)) + 
  geom_point(alpha=0.7)+ 
  scale_size(range = c(.1, 10)) +
  scale_colour_ipsum() +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("Taxa (Family)") +
  xlab("") +
  theme(legend.position = "none")

###pseq core ----

ggplot(pseq_core_larvae, aes(x=Treatment, y=OTU, size = Abundance, color = Treatment)) + 
  geom_point(alpha=0.7)+ 
  scale_size(range = c(.1, 10)) +
  scale_colour_ipsum() +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("Taxa (Family)") +
  xlab("") +
  theme(legend.position = "none")

###top phyla ----

ggplot(Top5P, aes(x=Age, y=OTU, size = Abundance, color = Age)) + 
  geom_point(alpha=0.7)+ 
  scale_size(range = c(.1, 10)) +
  scale_colour_ipsum() +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("Taxa (Family)") +
  xlab("") +
  theme(legend.position = "none")


#all samples

#core ASVs ----
#source: https://microbiome.github.io/tutorials/Core.html 

pseq <-  Marissa_MU42022_rare

###Day 1 only ----

pseq <- subset_samples(pseq, !Age %in% c("Spat", "Day 03", "Day 06", "Day 15"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

###psmelt ----

pseq <- psmelt(pseq)

pseq_core <- psmelt(pseq.core)

pseq_fam <- psmelt(pseq_fam)

pseq_phy <- psmelt(pseq_phy)

pseq_phy_rel <- psmelt(pseq.phy.rel)

pseq_fam_rel <- psmelt(pseq.fam.rel)

pseq_gen_rel <- psmelt(pseq.gen.rel)
