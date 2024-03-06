#Look at pseq.core

#Quick lib ----

library("devtools")
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


#Load data ----

Marissa_MU42022_rarefied_20231016 <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rarefied_20231016.rds")

pseq <-  Marissa_MU42022_rarefied_20231016

pseq <- Marissa_mb2021_filtered_20240203

#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

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


#remove F4 ----

pseq <- subset_samples(pseq, !Genetics %in% c("4"))

#Remove day 3 (only 1 sample remaining) for mb2021 project

pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))

#Day 1 only ----

pseq <- subset_samples(pseq, !Age %in% c("Spat", "Day 03", "Day 06", "Day 15"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#Spat only ----

pseq <- subset_samples(pseq, !Age %in% c("Day 01", "Day 03", "Day 06", "Day 15"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

pseq <- subset_samples(pseq, !Genetics %in% "4")

#remove algae ----

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#Larvae only with algae ----

pseq <- subset_samples(pseq, !Age %in% "Spat")

#Larvae only without algae ----

pseq <- subset_samples(pseq, !Age %in% "Spat")

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")


#Top phyla, all samples ----

top5F.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Family"], sum), TRUE)[1:5]

top5F = subset_taxa(pseq, Family %in% names(top5F.names))


top10G.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Genus"], sum), TRUE)[1:10]

top10G = subset_taxa(pseq, Genus %in% names(top10G.names))

#Create pseq objects ----


#convert to compositional data

pseq.rel <- microbiome::transform(pseq, "compositional")

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- microbiome::transform(pseq.core, "compositional")

pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq, level = "Genus", detection = 50/100, prevalence = 70/100)

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 95/100)


#pseq core - Probiotics + PB + Heat only (must remove all other samples first

pseq <- subset_samples(pseq, !Treatment %in% c("Control", "High temperature"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

pseq.core <- core(pseq, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")


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


#Top phyla psmelt ----

top5P <- psmelt(top5P)

top10G <- psmelt(top10G)

#If not using psmelt: 
#Extract core data as .csv ----


file_path <- "C:/Users/maris/OneDrive/Documents/GitHub/pseq_otu.csv"

# Write the otu_table to a CSV file
write.csv(pseq@otu_table, file = file_path)

write.csv(sam_data, file = file_path)

if (file.exists(file_path)) {
  cat("OTU table CSV file has been created:", file_path, "\n")
} else {
  cat("Failed to create the CSV file.")
}

#Have to manually change sample names into treatments in excel sheet
#import excel sheet - have to make it with x, y, z setup (data as rows)


####Box plots ----

#Want to look at abundance of ASVs 3, 88, 201 in spat data

selected_rows_ASVs <- subset(pseq, OTU %in% c("ASV201", "ASV88", "ASV3"))

View(selected_rows_ASVs)

#since ASV18 and 7 are 99.8% similar - combine into one OTU (most likely same bacteria)

combined_data <- selected_rows_ASVs %>%
  filter(OTU %in% c("ASV7", "ASV18")) %>%
  group_by(Age, Treatment) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )





#graph

p <- ggplot(combined_data, aes(x = Age, y = Avg_Abundance, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymin = Avg_Abundance - SD_Abundance, ymax = Avg_Abundance + SD_Abundance),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Probiotic Abundance - Spat", x = "Treatment", y = "Abundance") +
  theme(plot.title = element_text(hjust = 0.5))

print(p)

p <- ggplot(combined_data, aes(x = Age, y = Avg_Abundance, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymin = Avg_Abundance - SD_Abundance, ymax = Avg_Abundance + SD_Abundance),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Probiotic Abundance - Spat", x = "Treatment", y = "Abundance") +
  theme(plot.title = element_text(hjust = 0.5))

print(p)



p <- selected_rows_ASVs %>%
  mutate(name = fct_relevel(name, 
                            "ASV201", "ASV3", "ASV88")) %>%
  ggplot(selected_rows_ASVs, aes(x=OTU, y=Abundance, fill=Treatment)) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot() + scale_fill_brewer(palette = "Set2")
p


p <- ggplot(selected_rows_ASVs, aes(x=OTU, y=Abundance, fill=Treatment)) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot() + scale_fill_brewer(palette = "Set2" + xlab = "")
p

#with jitter - shows data distribution better than typical boxplot

p <- ggplot(selected_rows_ASVs, aes(x = Age, y = Abundance, fill = Treatment)) +
  geom_boxplot()
  
print(p)



#composition boxplots ----

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

###Bubbble plots ----


###pseq ----

ggplot(top10, aes(x=Treatment, y=Order, size = Abundance, color = Order)) + 
  geom_point(alpha=1)+ 
  scale_size(range = c(5, 15)) +
  scale_colour_brewer() +
  theme_ipsum() +
  facet_wrap(~Age)
  theme(legend.position="bottom") +
  ylab("Taxa (Family)") +
  xlab("") +
  theme(legend.position = "none")

###pseq relative ----

ggplot(pseq, aes(x = Treatment, y = Phylum, size = Abundance, color = Order)) + 
    geom_point(alpha = 0.75) + 
    scale_size(range = c(5, 15)) +
    scale_color_brewer(type = 'qual', palette = 'Paired') +  # Using the Set3 palette with 12 colors
    theme_classic() +
    facet_wrap(~Age) +
    theme(legend.position = "right") +
    ylab("Taxa (Order)") +
    xlab("") +
    theme(axis.text.x = element_text(size=11, angle=45, hjust=1))

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

###top families ----

  ggplot(top5F, aes(x = Treatment, y = Family, size = Abundance, color = Family)) + 
    geom_point(alpha = 0.75) + 
    scale_size(range = c(5, 15)) +
    scale_color_brewer(type = 'qual', palette = 'Paired') +  # Using the Set3 palette with 12 colors
    theme_classic() +
    facet_wrap(~Age) +
    theme(legend.position = "right") +
    ylab("Taxa (Order)") +
    xlab("") +
    theme(axis.text.x = element_text(size=11, angle=45, hjust=1))
  

#Top order ----
  
  ggplot(top10, aes(x = Treatment, y = Order, size = Abundance, color = Order)) + 
    geom_point(alpha = 0.75) + 
    scale_size(range = c(5, 15)) +
    scale_color_brewer(type = 'qual', palette = 'Paired') +  # Using the Set3 palette with 12 colors
    theme_classic() +
    facet_wrap(~Age) +
    theme(legend.position = "right") +
    ylab("Taxa (Order)") +
    xlab("") +
    theme(axis.text.x = element_text(size=11, angle=45, hjust=1))
  
  #Top class ----
  
  ggplot(top10C, aes(x = Treatment, y = Class, size = Abundance, color = Class)) + 
    geom_point(alpha = 0.75) + 
    scale_size(range = c(5, 15)) +
    scale_color_brewer(type = 'qual', palette = 'Paired') +  # Using the Set3 palette with 12 colors
    theme_classic() +
    facet_wrap(~Age) +
    theme(legend.position = "right") +
    ylab("Taxa (Class)") +
    xlab("") +
    theme(axis.text.x = element_text(size=11, angle=45, hjust=1))  
  
  #Top order ----
  
  ggplot(top10, aes(x = Treatment, y = Order, size = Abundance, color = Order)) + 
    geom_point(alpha = 0.75) + 
    scale_size(range = c(5, 15)) +
    scale_color_brewer(type = 'qual', palette = 'Paired') +  # Using the Set3 palette with 12 colors
    theme_classic() +
    facet_wrap(~Age) +
    theme(legend.position = "right") +
    ylab("Taxa (Order)") +
    xlab("") +
    theme(axis.text.x = element_text(size=11, angle=45, hjust=1))    
  
#View Vibrionaceae abundance ----

selected_rows_ASVs <- subset(pseq, Family %in% c("Vibrionaceae"))

View(selected_rows_ASVs)

#combine all vibrio abundances

combined_data <- selected_rows_ASVs %>%
  filter(Family %in% c("Vibrionaceae")) %>%
  group_by(Age, Treatment) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )

View(combined_data)

#graph

p <- ggplot(combined_data, aes(x = Age, y = Avg_Abundance, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymin = Avg_Abundance - SD_Abundance, ymax = Avg_Abundance + SD_Abundance),
    position = position_dodge(width = 0.9),
    width = 0.25) +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Vibrionaceae Abundance", x = "", y = "Abundance") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(combined_data, aes(x = Age, y = Avg_Abundance, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Vibrionaceae", x = "", y = "Relative Abundance (Reads)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid = element_blank())


View(data)


#View Rhodo abundance ----

selected_rows_ASVs <- subset(pseq, Family %in% c("Rhodobacteraceae"))

View(selected_rows_ASVs)


p <- ggplot(selected_rows_ASVs, aes(x = Treatment, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") + facet_grid(cols = vars(Age))

print(p)

#make averages of genus's for each treatment and Age point

combined_data <- selected_rows_ASVs %>%
  group_by(Genus, Treatment, Age) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )

View(combined_data)

#graph

p <- ggplot(combined_data, aes(x = Treatment, y = Avg_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") + facet_grid(cols = vars(Age)) +
  geom_errorbar(
    aes(ymin = Avg_Abundance - SD_Abundance, ymax = Avg_Abundance + SD_Abundance),
    position = position_dodge(width = 0.9),
    width = 0.25) +
  labs(title = "Rhodobacteraceae Abundance", x = "", y = "Abundance") +
  theme(plot.title = element_text(hjust = 0.5))

print(p)

combined_data_2 <- subset(combined_data, Genus %in% c("Celeribacter"))

ggplot(combined_data_2, aes(x = Treatment, y = Avg_Abundance, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~Age) +
  geom_errorbar(
    aes(ymin = Avg_Abundance - SD_Abundance, ymax = Avg_Abundance + SD_Abundance),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +  scale_fill_brewer(palette = "Set2") + 
  labs(title = "Celeribacter Abundance", x = "", y = "Abundance") +
  theme(plot.title = element_text(hjust = 0.5))

