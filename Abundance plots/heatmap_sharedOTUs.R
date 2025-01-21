#Heat map of top ASVs
#trying to replicate heatmap done in Domin et al., 2023
#source: https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_extract_otus.R

devtools::install_github("vmikk/metagMisc")
install.packages("metagMisc")

library("devtools")
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)
library(metagMisc)
library(pheatmap)
library(tidyr)

#set theme
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

#Select ASVs "which were constantly detected in at least one of the four 
#different timepoints of recolonization across all samples (minimum relative abundance was 0.005%)" (Domin et al., 2023)

pseq <-`Filtered_Rarified_MU42022_23-12-13`
pseq <- MU42022_filtered_Oct92024
pseq <- PB2023_spat_not_rarefied_normalized

#Exclude factors

pseq <- subset_samples(pseq, Age %in% c("Spat"))

pseq <- subset_samples(pseq, !Treatment %in% c("James", "Continuous Probiotics"))

pseq <- subset_samples(pseq, !Organism %in% c("Algae"))


#check if any OTUs are not present in any samples (want false)
any(taxa_sums(pseq) == 0)

#if true
pseq_filtered <- prune_taxa(taxa_sums(pseq) > 0, pseq)
any(taxa_sums(pseq_filtered) == 0)
pseq <- pseq_filtered

# Transform to relative abundance
pseq_rel <- transform_sample_counts(pseq, function(x) x / sum(x))
# Subset to top 20 ASVs
# top_asvs <- names(sort(taxa_sums(pseq_rel), decreasing = TRUE))[1:579]
top_asvs <- names(sort(taxa_sums(pseq_rel), decreasing = TRUE))[1:141]
pseq_top <- prune_taxa(top_asvs, pseq_rel)

# Generate a heatmap
plot_heatmap(pseq_top, 
             sample.label = "Treatment", 
             taxa.label = F, 
             low = "white", 
             high = "red", 
             na.value = "gray")


#Taxa shared between all samples samples

#source: https://microbiome.github.io/tutorials/Core.html

#Merge ASv data for each time point
#Convert to relative comp
pseq_filt <- microbiome::transform(pseq_filt, "compositional")

ps10<-merge_samples(pseq_filt, "Age", fun= mean)
ps10@otu_table[,1:10]
?core_members
#Select core
Core_trial<-core(ps10, detection=0.5/100, prevalence = 90/100)
View(Core_trial@otu_table)
Core_trial


Core_trial <- psmelt(Core_trial)
str(Core_trial)

#create list of ASV names

otu_names_list <- Core_trial$OTU
otu_names_list

#recall pseq_filt object
pseq_filt

#prune taxa to only include taxa on list (shared between treatments)

otu_names_list <- Core_trial$OTU

prune_pseq_filt<- prune_taxa(otu_names_list, pseq_filt)

prune_pseq_filt

prune_psmelt <- psmelt(prune_pseq_filt)


#All data (not just shared otus)
pseq <- microbiome::transform(pseq, "compositional")
prune_psmelt <- psmelt(pseq)
  


# Order Age as a factor in the desired order
# prune_psmelt$Age <- factor(prune_psmelt$Age, levels = c("Day 01", "Day 03", "Day 06", "Day 15", "Spat"))
#prune_psmelt$Age <- factor(prune_psmelt$Day, levels = c("Day 01", "Day 03", "Day 06", "Day 15", "Spat"))


# Summarize abundance by OTU and Age
heatmap_data <- prune_psmelt %>%
  group_by(OTU, Treatment) %>%
  summarize(Mean_Abundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = Treatment, values_from = Mean_Abundance, values_fill = 0)

# Convert the OTU column into row names for clustering
heatmap_matrix <- as.data.frame(heatmap_data)
rownames(heatmap_matrix) <- heatmap_matrix$OTU
heatmap_matrix$OTU <- NULL  # Remove OTU column since it's now row names

# Scale the data by row (optional: for clearer visualization)
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))  # Row-wise scaling

#If getting errors
anyNA(heatmap_matrix_scaled) # Check for NA values
any(is.nan(heatmap_matrix_scaled)) # Check for NaN values
any(is.infinite(heatmap_matrix_scaled)) # Check for infinite values
# Replace NA and NaN with 0
heatmap_matrix_scaled[is.na(heatmap_matrix_scaled)] <- 0
heatmap_matrix_scaled[is.nan(heatmap_matrix_scaled)] <- 0



# Plot the heatmap

pheatmap(heatmap_matrix_scaled,
         color = colorRampPalette(c("skyblue2", "white", "red"))(50), # Custom gradient
         cluster_rows = TRUE,   # Clustering OTUs
         cluster_cols = FALSE,  # No clustering for Age
         show_rownames = FALSE, # Optional: hide OTU names
         fontsize_col = 18,
         angle_col = 0,
         main = "")



#pheatmaps with clusters ----
library(dplyr)
library(tidyr)
library(pheatmap)

# Prepare data (from your previous steps)
heatmap_data <- prune_psmelt %>%
  group_by(OTU, Treatment) %>%
  summarize(Mean_Abundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = Treatment, values_from = Mean_Abundance, values_fill = 0)

# Convert to matrix for pheatmap
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) <- heatmap_data$OTU

# Scale the data
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))  # Row-wise scaling

#If getting errors
anyNA(heatmap_matrix_scaled) # Check for NA values
any(is.nan(heatmap_matrix_scaled)) # Check for NaN values
any(is.infinite(heatmap_matrix_scaled)) # Check for infinite values
# Replace NA and NaN with 0
heatmap_matrix_scaled[is.na(heatmap_matrix_scaled)] <- 0
heatmap_matrix_scaled[is.nan(heatmap_matrix_scaled)] <- 0

# Perform hierarchical clustering
row_hclust <- hclust(dist(heatmap_matrix_scaled), method = "complete")

# Cut the tree into clusters (e.g., 5 clusters)
k <- 9  # Number of clusters
row_clusters <- cutree(row_hclust, k = k)

# Create a data frame for row annotations
row_annotation <- data.frame(Cluster = factor(row_clusters))
rownames(row_annotation) <- rownames(heatmap_matrix_scaled)

# Assign colors for the clusters
annotation_colors <- list(Cluster = c(
  "1" = "yellow3",
  "2" = "orange",
  "3" = "forestgreen",
  "4" = "red4",
  "5" = "pink",
  "6" = "grey",
  "7" = "red3",
  "8" = "green",
  "9" = "purple",
  "10" = "pink3",
  "11" = "blue",
  "12" = "limegreen"
))


# Generate the heatmap with cluster labels
pheatmap(
  heatmap_matrix_scaled,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,  # Ensure clustering
  cluster_cols = FALSE,  # No clustering on columns
  show_rownames = FALSE,  # Optionally hide row names for clarity
  main = "",
  angle_col = 0,
  fontsize_col = 18,
  fontsize = 14
)

#Tree of ASV clusters ----
library(ape)
library(ggtree)
library(ggplot2)

# Convert hclust to phylo object
tree <- as.phylo(row_hclust)
# Create a data frame with ASV and cluster information
tip_data <- data.frame(
  label = names(row_clusters),  # ASV names
  Cluster = factor(row_clusters)  # Cluster assignment
)

# Define cluster colors
cluster_colors <- annotation_colors$Cluster

# Create the circular tree
 ggtree(tree, layout = "fan", branch.length = 0.5) %<+% tip_data +
  geom_tiplab(aes(color = Cluster), size = 5, show.legend = FALSE) +
  geom_point(aes(color = Cluster), size = 1.5) +  # Add points to tips
  scale_color_manual(values = cluster_colors) +
  theme(legend.position = "right")


#Look at cluster 2 only 
# Filter rows for cluster 1
cluster_1_rows <- rownames(row_annotation)[row_annotation$Cluster == "3"]
heatmap_matrix_cluster_1 <- heatmap_matrix_scaled[cluster_1_rows, , drop = FALSE]

# Subset the row annotation to match the filtered rows
row_annotation_cluster_1 <- row_annotation[cluster_1_rows, , drop = FALSE]

# Generate the heatmap for cluster 1
pheatmap(
  heatmap_matrix_cluster_1,
  annotation_row = row_annotation_cluster_1,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,  # Ensure clustering within cluster 1 rows
  cluster_cols = FALSE,  # No clustering on columns
  show_rownames = TRUE,  # Show row names for cluster 1
  main = "Cluster 3 Heatmap",
  angle_col = 0,
  fontsize_col = 18,
  fontsize_row = 14,
  fontsize = 14
)


#Look at ASVs in each cluster ----
#Cluster 1 looks ifferent in spat between control/PB and PBH

pseq <- MU42022_filtered_Oct92024
# pseq_filt <- subset_samples(pseq, Treatment %in% c("Control"))
pseq_filt <- subset_samples(pseq_filt, Age %in% c("Spat"))
#pseq_filt <- microbiome::transform(pseq_filt, "compositional")
prune_psmelt <- psmelt(pseq_filt)

# Create a data frame with ASVs and their corresponding cluster assignments
cluster_summary <- data.frame(ASV = rownames(heatmap_matrix_scaled),
                              Cluster = row_clusters)

# Summarize ASVs by cluster
cluster_1 <- cluster_summary %>%
  group_by(Cluster) %>%
  filter(Cluster %in% c("1"))

# View the cluster details
#View(cluster_details)

# Extract the list of ASVs from cluster_1
cluster_1_asvs <- cluster_1$ASV

# Filter prune_psmelt to include only rows with ASVs in cluster_1
prune_psmelt_cluster1 <- prune_psmelt %>% 
  filter(OTU %in% cluster_1_asvs)

Avg_abundance <- prune_psmelt_cluster1 %>%
  group_by(Treatment, Order, Family, Genus) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  ) %>%
  group_by(Treatment) %>%
  mutate(Avg_Abundance = 100 * Avg_Abundance / sum(Avg_Abundance))

paired_palette <- brewer.pal(12, "Paired")

extended_palette <- c(paired_palette, "pink",  "yellow", "lightgreen", "green", "brown", "orange", "red")  # Add a custom color to the palette

# Use ggplot with the extended Paired palette
ggplot(Avg_abundance, aes(fill = Order, y = Avg_Abundance, x = Treatment)) + 
  geom_bar(position = "stack", stat = "identity", colour = "black") +
  scale_fill_manual(values = extended_palette) +  # Use scale_fill_manual to specify the extended palette
  labs(title = "", x = "", y = "Relative abundance (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1, face = "bold"),
    axis.title.y = element_text(size = 14),
    strip.background = element_rect(fill = "white"),  # Set facet heading background to white
    strip.text = element_text( size = 12)           # Set facet heading text to bold
  )





# Step 1: Prepare the data
heatmap_data <- prune_psmelt %>%
  group_by(OTU, Sample.ID) %>% 
  summarise(Abundance = sum(Abundance), .groups = "drop") %>% # Summarize if necessary
  pivot_wider(names_from = Sample.ID, values_from = Abundance, values_fill = 0) # Wide format

# Convert to a matrix
heatmap_matrix <- as.data.frame(heatmap_data)
rownames(heatmap_matrix) <- heatmap_matrix$OTU
heatmap_matrix <- heatmap_matrix[, -1] # Remove OTU column

# Scale the data by row (optional: for clearer visualization)
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))  # Row-wise scaling

#If getting errors
anyNA(heatmap_matrix_scaled) # Check for NA values
any(is.nan(heatmap_matrix_scaled)) # Check for NaN values
any(is.infinite(heatmap_matrix_scaled)) # Check for infinite values
# Replace NA and NaN with 0
heatmap_matrix_scaled[is.na(heatmap_matrix_scaled)] <- 0
heatmap_matrix_scaled[is.nan(heatmap_matrix_scaled)] <- 0

# Generate the heatmap
pheatmap(heatmap_matrix_scaled,
         color = colorRampPalette(c("skyblue2", "white", "red"))(50), # Custom gradient
         cluster_rows = TRUE,  
         cluster_cols = FALSE,  
         show_rownames = FALSE, 
         fontsize_col = 12,    
         main = "Control")

#replace x-axis with treatment names
# Extract the metadata
metadata <- as.data.frame(pseq@sam_data)

# Define the desired treatment order
treatment_order <- c("Control", "Killed-Probiotics", "Probiotics")

# Sort metadata by treatment order
metadata <- metadata[order(factor(metadata$Treatment, levels = treatment_order)), ]

# Reorder the columns of the heatmap matrix to match the sorted metadata
heatmap_matrix_ordered <- heatmap_matrix_scaled[, metadata$Sample.ID]

# Use the Treatment variable as custom x-axis labels
x_labels <- metadata$Treatment

# Generate the heatmap with ordered samples and custom labels
pheatmap(heatmap_matrix_ordered,
         color = colorRampPalette(c("skyblue2", "white", "red"))(50), # Custom gradient
         cluster_rows = TRUE,   # Clustering OTUs
         cluster_cols = FALSE,  # No clustering for columns
         show_rownames = FALSE, # Optional: hide OTU names
         main = "Heatmap by Treatment",
         angle_col = 45,
         fontsize_col = 14,
         fontsize_row = 14,
         fontsize = 14,
         labels_col = x_labels)


