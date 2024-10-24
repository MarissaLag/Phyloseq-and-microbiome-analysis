#plot ASVs overtime

library(phyloseq)
library(ggplot2)
library(dplyr)
library(microbiome)

#MU42022 filtering
pseq <- MU42022_filtered_Oct92024
#pseq <- MU42022_filtered_NOT_rarefied
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")

# List of ASVs you're interested in
asv_ids <- c("ASV88", "ASV201")
asv_ids <- c("ASV7", "ASV18")
asv_ids <- c("ASV11", "ASV198", "ASV201", "ASV471", "ASV613")
asv_ids <- c("ASV444")

# Filter the phyloseq object to only include these ASVs
pseq_filtered <- prune_taxa(taxa_names(pseq) %in% asv_ids, pseq)
pseq_filtered <- psmelt(pseq_filtered)

pseq_filtered$OTU <- factor(pseq_filtered$OTU, levels = c("ASV88", "ASV201"))
pseq_filtered$OTU <- factor(pseq_filtered$OTU, levels = c("ASV7", "ASV18"))
pseq_filtered$OTU <- factor(pseq_filtered$OTU, levels = c("ASV11", "ASV198", "ASV201", "ASV471", "ASV613"))

ggplot(pseq_filtered, aes(x = Age, y = Abundance, color = Treatment)) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(OTU, Treatment)), linetype = "dashed") +  # Dashed lines to connect points
  facet_wrap(~ OTU, scales = "free_y") +  # Facet by OTU with custom ordering
  scale_color_manual(values = c("darkgrey", "cornflowerblue", "orange")) +  # Custom color for Treatment
  labs(title = "",
       x = "Age",
       y = "Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank())

#Average abundace for each treatmentxage group
pseq_avg <- pseq_filtered %>%
  group_by(OTU, Treatment, Age, Genus) %>%
  summarise(Avg_Abundance = mean(Abundance), .groups = 'drop')

ggplot(pseq_avg, aes(x = Age, y = Avg_Abundance, color = Treatment)) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(OTU, Treatment)), linetype = "dashed") +  # Dashed lines to connect points
  facet_wrap(~ otu.rare, scales = "free_y") +  # Facet by OTU with custom ordering
  scale_color_manual(values = c("darkgrey", "cornflowerblue", "orange")) +  # Custom color for Treatment
  labs(title = "",
       x = "",
       y = "Average Abundance") +  # Change Y-axis label to indicate average
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Larger and slanted x-axis labels
        strip.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size = 14))  # Larger facet titles

