#Abundance plots using psmelt
#Written by MArissa WL (2024-01-18)
#workflow: Load data (pseq) -> filter data (if you want) -> manipulate data (e.g., top families, relative composition, etc.) -> psmelt -> plot


library("devtools")
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)
library(dplyr)

#Load data ----

Marissa_MU42022_rarefied_20231016 <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rarefied_20231016.rds")

pseq <-  Marissa_MU42022_rarefied_20231016

#set theme ----

theme.marissa <- function() {
  theme_classic(base_size = 14) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"))
}

theme_set(theme.marissa())

#Sample selection ----

#Remove F4 ----

pseq <- subset_samples(pseq, !Genetics %in% c("4"))

#Remove algae ----

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#Select time-points ----

pseq <- subset_samples(pseq, !Age %in% c("Spat", "Day 03", "Day 06", "Day 15"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#Spat only

pseq <- subset_samples(pseq, !Age %in% c("Day 01", "Day 03", "Day 06", "Day 15"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

pseq <- subset_samples(pseq, !Genetics %in% "4")


#Larvae only with algae

pseq <- subset_samples(pseq, !Age %in% "Spat")


#Top families ----

top5F.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Family"], sum), TRUE)[1:5]

top5F = subset_taxa(pseq, Family %in% names(top5F.names))

#Top Genus ----

top10G.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Genus"], sum), TRUE)[1:10]

top10G = subset_taxa(pseq, Genus %in% names(top10G.names))

#Relative composition ----

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- microbiome::transform(pseq.core, "compositional")

pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq, level = "Genus", detection = 50/100, prevalence = 70/100)

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 95/100)


#psmelt ----

#all data

top5F <- psmelt(top5F)

top10G <- psmelt(top10G)

pseq <- psmelt(pseq)


#composition boxplots ----

Avg_abundance <- top5F %>%
  group_by(Treatment, Family) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )

Avg_abundance <- pseq %>%
  group_by(Treatment, Family) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )

ggplot(Avg_abundance, aes(fill=Family, y=Avg_Abundance, x=Treatment)) + 
  geom_bar(position="dodge", stat="identity") + scale_fill_ipsum()

ggplot(Avg_abundance, aes(fill=Family, y=Avg_Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_brewer(palette = "Spectral") +

ggplot(top5F, aes(fill=Family, y=Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Spat samples", x = "", y = "Average abundance") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(Avg_abundance, aes(fill=Family, y=Avg_Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_ipsum() +
  labs(title = "Larval samples", x = "", y = "Average abundance") +
  theme(plot.title = element_text(hjust = 0.5)) 

#Genus level ----
Avg_abundance <- top10G %>%
  group_by(Treatment, Genus) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )

ggplot(Avg_abundance, aes(fill=Genus, y=Avg_Abundance, x=Treatment)) + 
  geom_bar(position="dodge", stat="identity")

ggplot(Avg_abundance, aes(fill=Genus, y=Avg_Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity")
