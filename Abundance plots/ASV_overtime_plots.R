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
pseq <- microbiome::transform(pseq, "compositional")


# List of ASVs you're interested in
asv_ids <- c("ASV88", "ASV178")
asv_ids <- c("ASV7", "ASV18")
asv_ids <- c("ASV11", "ASV198", "ASV201", "ASV471", "ASV613")
asv_ids <- c("ASV316")

# Filter the phyloseq object to only include these ASVs
pseq_filtered <- prune_taxa(taxa_names(pseq) %in% asv_ids, pseq)
pseq_filtered <- psmelt(pseq_filtered)

pseq_filtered$OTU <- factor(pseq_filtered$OTU, levels = c("ASV88", "ASV178"))
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
  group_by(OTU, Treatment, Age) %>%
  summarise(Avg_Abundance = mean(Abundance), .groups = 'drop')


# Create a named vector to map old facet labels to new ones
# facet_labels <- c(
#   "ASV88" = "ASV88 - Loktanella",
#   "ASV178" = "ASV178 - Litoreibacter"
# )

facet_labels <- c(
  "ASV316" = "ASV316 - Colwelliaceae"
)

# facet_labels <- c(
#   "ASV11" = "ASV11 - Uncultured Roseobacter",
#   "ASV198" = "ASV198 - Uncultured Roseobacter",
#   "ASV201" = "ASV201 - Sulfitobacter",
#   "ASV471" = "ASV471 - Sulfitobacter",
#   "ASV613" = "ASV613 - Uncultured Roseobacter"
# )

# Use the labeller argument in facet_wrap()
ggplot(pseq_avg, aes(x = Age, y = Avg_Abundance, color = Treatment)) +
  geom_point(size = 4) +
  geom_line(aes(group = interaction(OTU, Treatment)), linetype = "dashed") +  
  facet_wrap(~ OTU, scales = "free_y", labeller = labeller(OTU = facet_labels)) +  
  scale_color_manual(values = c("darkgrey", "cornflowerblue", "orange")) +  
  labs(title = "",
       x = "",
       y = "Average Relative Abundance") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white", color = "black"), # White background for facet headings
    strip.text = element_text(size = 14, face = "bold"), # Customize text
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

#combined OTUs
ggplot(pseq_avg, aes(x = Age, y = Avg_Abundance, color = Treatment)) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(Treatment)), linetype = "dashed") +  # Dashed lines to connect points
  #facet_wrap(~ otu.rare, scales = "free_y") +  # Facet by OTU with custom ordering
  scale_color_manual(values = c("darkgrey", "cornflowerblue", "orange")) +  # Custom color for Treatment
  labs(title = "",
       x = "",
       y = "Average Relative Abundance") +  # Change Y-axis label to indicate average
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Larger and slanted x-axis labels
        strip.background = element_rect(fill = "white", color = "black"),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size = 14)) 
