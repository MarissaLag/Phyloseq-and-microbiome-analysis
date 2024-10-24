#Microbiome composition plots - correcting for uneven sample numbers between groups

library("devtools")
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

#MU42022 filtering
pseq <- MU42022_filtered_Oct92024
#pseq <- MU42022_filtered_NOT_rarefied
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")
pseq <- subset_samples(pseq, Age %in% c("Spat"))

#MB2021 filtering
pseq <- Marissa_mb2021_filtered_20240203
pseq <- mb2021_filtered_NOT_rarefied_normalized
pseq <- subset_samples(pseq, Age %in% c("Day 01"))
pseq <- subset_samples(pseq, !Family %in% c("9")) #remove T9 spat samples


#1st create new column that contains Treat*Age information 
sample_data <- pseq@sam_data
# Create a new column combining Age and Treatment ----
sample_data$Age_Treatment <- paste(sample_data$Age, sample_data$Treatment, sep = "_")
sample_data(pseq) <- sample_data
View(pseq@sam_data)

#Convert to relative abundance ----
#result from either transformation the same after averaging

#percentage (out of 100%)
pseq2 = filter_taxa(pseq, function(x) mean(x) > 0.1, TRUE)
pseq2

pseq3 = transform_sample_counts(pseq2, function(x) x / sum(x) )
pseq3

#other relative transformation
pseq3 <- microbiome::transform(pseq, "compositional")

#Top 10 families

top10F.names = sort(tapply(taxa_sums(pseq3), tax_table(pseq3)[, "Family"], sum), TRUE)[1:10]

top10F = subset_taxa(pseq3, Family %in% names(top10F.names))

#Top 10 genera 
top10G.names = sort(tapply(taxa_sums(pseq3), tax_table(pseq3)[, "Genus"], sum), TRUE)[1:10]

top10G = subset_taxa(pseq3, Genus %in% names(top10G.names))

#core
pseq.core <- core(pseq3, detection = .1/100, prevalence = 90/100)

#rare
pseq <- aggregate_rare(pseq3, level = "Genus", detection = 1/100, prevalence = 50/100)

#top genera of core
top10G.names = sort(tapply(taxa_sums(pseq.core), tax_table(pseq.core)[, "Genus"], sum), TRUE)[1:10]

top10G.core = subset_taxa(pseq.core, Genus %in% names(top10G.names))


#psmelt ----
#Convert to data frame
pseq_psmelt <- psmelt(top10F)
pseq_psmelt <- psmelt(top10G)
pseq_psmelt <- psmelt(pseq)
pseq_psmelt <- psmelt(pseq.core)
pseq_psmelt <- psmelt(top10G.core)


#plot of rel abundaces of each sample (out of 100%)
ggplot(pseq_psmelt, aes(fill=Phylum, y=Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_brewer(palette = "Paired") +
  labs(title = "All time-points", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank())

#look at sample numbers per group ----

sample_counts <- pseq_psmelt %>%
  group_by(Treatment) %>%
  summarise(Num_Samples = n())

print(sample_counts)

#results of sample numbers
#each sample has Num_Samples of 580 (i.e., number of taxa)
#Day 1 samples - Control/HS/LS = 3
#Day 18 samples - Control/HS = 4, LS = 3
#Spat samples - Control = 10, HS = 12, LS = 10


#Average abundance per sample group ----

Avg_abundance <- pseq_psmelt %>%
  group_by(Age_Treatment, Family) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  ) %>%
  separate(Age_Treatment, into = c("Age", "Treatment"), sep = "_") 

View(Avg_abundance)

#Make abundance out of 100%?
Avg_abundance <- pseq_psmelt %>%
  group_by(Age_Treatment, Genus) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  ) %>%
  group_by(Age_Treatment) %>%
  mutate(Avg_Abundance = 100 * Avg_Abundance / sum(Avg_Abundance)) %>% 
separate(Age_Treatment, into = c("Age", "Treatment"), sep = "_") 





#plot ----

paired_palette <- brewer.pal(12, "Paired")

# Add another color to the Paired palette
extended_palette <- c(paired_palette, "#ff7f00", "pink", "red", "yellow", "lightgreen")  # Add a custom color to the palette

# Use ggplot with the extended Paired palette
p3 <- ggplot(Avg_abundance, aes(fill = Genus, y = Avg_Abundance, x = Treatment)) + 
  geom_bar(position = "stack", stat = "identity", colour = "black") +
  scale_fill_manual(values = extended_palette) +  # Use scale_fill_manual to specify the extended palette
  labs(title = "", x = "", y = "Relative abundance (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1, face = "bold")) +
  facet_wrap("Age") +
  theme(legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.title.x = element_text(size = 5),
        axis.title.y = element_text(size = 14))
p3  


#Vibrionaceae abundances ----
pseq <- mb2021_filtered_NOT_rarefied
pseq<- Marissa_mb2021_filtered_20240203
pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))
pseq <- subset_samples(pseq, !Family %in% c("9"))

pseq <- subset_samples(pseq, Age %in% c("Spat")) 

#correct family column

pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

View(pseq@sam_data)

#1st create new column that contains Treat*Age*Family information 
sample_data <- pseq@sam_data
# Create a new column combining Age and Treatment ----
sample_data$Age_Treatment_Family <- paste(sample_data$Age, sample_data$Treatment, sample_data$Family, sep = "_")
sample_data(pseq) <- sample_data
View(pseq@sam_data)

pseq3 <- microbiome::transform(pseq, "compositional")

#remove vibrio after or before compositional transformation

pseq3 <- subset_taxa(pseq3, Family=="Vibrionaceae")

pseq_psmelt <- psmelt(pseq3)

View(pseq_psmelt)

Avg_abundance <- pseq_psmelt %>%
  group_by(Age_Treatment_Family, Genus) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  ) %>%
  group_by(Age_Treatment_Family) %>%
  mutate(Avg_Abundance = 100 * Avg_Abundance / sum(Avg_Abundance)) %>% 
  separate(Age_Treatment_Family, into = c("Age", "Treatment"), sep = "_") 

Avg_abundance <- pseq_psmelt %>%
  group_by(Age_Treatment_Family, Genus) %>%
  summarise(
    Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  ) %>%
  separate(Age_Treatment_Family, into = c("Age", "Treatment", "Family"), sep = "_") 


sum_abundance <- pseq_psmelt %>%
  group_by(Age_Treatment, Genus) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = 'drop'
  ) %>% 
  separate(Age_Treatment, into = c("Age", "Treatment"), sep = "_")

ggplot(Avg_abundance, aes(fill = Genus, y = Abundance, x = Treatment)) + 
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = extended_palette) +
  labs(title = "Vibrionaceae (T9 removed)", 
       x = "", 
       y = "Averaage Relative Abundance", 
       fill = "Vibrionaceae (Genus)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),  # Adjust y-axis title font size
        legend.text = element_text(size = 12),   # Adjust legend text font size
        legend.title = element_text(size = 14),
        panel.border = element_blank()) +  # Adjust legend title font size
  facet_wrap("Age")


ggplot(Avg_abundance, aes(fill = Genus, y = Abundance, x = Treatment)) + 
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = extended_palette) +
  labs(title = "Vibrionaceae (T9 removed)", 
       x = "", 
       y = "Averaage Relative Abundance", 
       fill = "Vibrionaceae (Genus)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),  # Adjust y-axis title font size
        legend.text = element_text(size = 12),   # Adjust legend text font size
        legend.title = element_text(size = 14),
        panel.border = element_blank()) +  # Remove border lines around the panels
  facet_grid(Age ~ Family)


  