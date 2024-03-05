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

#pseq<- Marissa_mb2021_filtered_20240203
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

#Remove F4 (MU42022) ----

pseq <- subset_samples(pseq, !Genetics %in% c("4"))

#Remove day 3 (only 1 sample remaining) for mb2021 project

pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))

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

pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq, level = "Genus", detection = 50/100, prevalence = 70/100)

pseq.core <- core(pseq_fam, detection = .1/100, prevalence = 95/100)


pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- microbiome::transform(pseq.core, "compositional")

#Top families relative abundance ----

top5F.names = sort(tapply(taxa_sums(pseq.fam.rel), tax_table(pseq.fam.rel)[, "Family"], sum), TRUE)[1:10]

top5F = subset_taxa(pseq.fam.rel, Family %in% names(top5F.names))

#Top Genus relative abundance ----

top10G.names = sort(tapply(taxa_sums(pseq.gen.rel), tax_table(pseq.gen.rel)[, "Genus"], sum), TRUE)[1:10]

top10G = subset_taxa(pseq.gen.rel, Genus %in% names(top10G.names))

#psmelt ----

#all data

top5F <- psmelt(top5F) 

top10G <- psmelt(top10G)

pseq <- psmelt(pseq)

relative_fam <- psmelt(pseq.fam.rel)

relative_fam <- psmelt(pseq.phy.rel)

core <- psmelt(pseq.core)

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
  labs(title = "All time-points", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(top10G, aes(fill=Genus, y=Abundance, x=Treatment)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "All time-points", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(relative, aes(fill=Family, y=Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "All time-points", x = "", y = "Total abundance") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(core, aes(fill=Family, y=Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Core - All time-points", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Age) +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
  theme(axis.text.y = element_text(size=11)) +
  theme(legend.text = element_text(size = 11))
  
  

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
