#Phylogenetic relatedness between ASVs
library("devtools")
library(phyloseq)
library(microbiome)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(ggtree)
library(ape)

MU42022_filtered_Oct92024 <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/MU42022_filtered_Oct92024.rds")

pseq <- MU42022_filtered_Oct92024
pseq <- microbiome::transform(pseq, "compositional")
pseq <- psmelt(pseq) 

pseq <- PB2023_spat_not_rarefied_normalized
pseq_filt <- subset_samples(pseq, !Treatment %in% c("Continuous Probiotics", "James"))
any(taxa_sums(pseq_filt) == 0)
# pseq_filtered <- prune_taxa(taxa_sums(pseq_filt) > 0, pseq_filt)
# any(taxa_sums(pseq_filtered) == 0)
pseq_filt <- pseq_filtered
pseq <- microbiome::transform(pseq_filt, "compositional")
pseq <- psmelt(pseq)


#Merge ASV seq with metadata
colnames(MU42022_sequence_ASVname_mapping) <- c("ASV", "Sequence")
# colnames(sequence_ASVname_mapping_SMK) <- c("ASV", "Sequence")
# head(sequence_ASVname_mapping_SMK)

# Merge the ASV table with the sequence data based on ASV number
merged_data <- merge(pseq, MU42022_sequence_ASVname_mapping, by.x = "OTU", by.y = "ASV", all.x = TRUE)
# merged_data <- merge(pseq, sequence_ASVname_mapping_SMK, by.x = "OTU", by.y = "ASV", all.x = TRUE)


#Trees ----
#Subset data for Rhodobacteraceae family
rhodobacteraceae_data <- subset(merged_data, Family == "Rhodobacteraceae")
#Extract unique sequences for the Rhodobacteraceae ASVs
sequences <- unique(rhodobacteraceae_data$Sequence)
names(sequences) <- unique(rhodobacteraceae_data$OTU)  # Assign ASV names to the sequences

# Subset data for alphaproteobacteria 
alpha_data <- subset(merged_data, Class == "Alphaproteobacteria")
#View(alpha_data)
# Extract unique sequences for the Rhodobacteraceae ASVs
sequences <- unique(alpha_data$Sequence)
names(sequences) <- unique(alpha_data$OTU)  # Assign ASV names to the sequences

#All data 
# Extract unique sequences for the Rhodobacteraceae ASVs
sequences <- unique(merged_data$Sequence)
names(sequences) <- unique(merged_data$OTU)

dna_sequences <- DNAStringSet(sequences)

# Align the sequences
aligned_sequences <- AlignSeqs(dna_sequences)
# Convert DNAStringSet to a character matrix
aligned_matrix <- as.matrix(aligned_sequences)

# Convert the matrix to a phyDat object
aligned_phyDat <- phyDat(aligned_matrix, type = "DNA")

# Check the phyDat object
print(aligned_phyDat)

# Compute the distance matrix
dist_matrix <- dist.ml(aligned_phyDat)

# Build the phylogenetic tree
tree <- NJ(dist_matrix)

# Plot the tree
# Plot the tree with smaller ASV labels
plot(tree, main = "Phylogenetic Tree", cex = 0.35)  # Adjust cex as needed

highlighted_ASVs_red <- c("ASV11", "ASV88", "ASV198", "ASV178", "ASV201", "ASV471", "ASV613")
highlighted_ASVs_blue <- c("ASV7", "ASV18")

#PB2023
highlighted_ASVs_red <- c("ASV190")
highlighted_ASVs_blue <- c("ASV227")

# Get the tip labels
tip_labels <- tree$tip.label

# Add red and blue bold labels for the specified ASVs
for (i in seq_along(tip_labels)) {
  if (tip_labels[i] %in% highlighted_ASVs_red) {
    tiplabels(text = tip_labels[i], tip = i, col = "red", font = 2, cex = 1)
  } else if (tip_labels[i] %in% highlighted_ASVs_blue) {
    tiplabels(text = tip_labels[i], tip = i, col = "blue", font = 2, cex = 1)
  }
}

#remove all other tip labels
# Get the tip labels
tip_labels <- tree$tip.label

# Create a vector for labels, only keeping the selected ASVs and removing the rest
new_tip_labels <- ifelse(tip_labels %in% c(highlighted_ASVs_red, highlighted_ASVs_blue), tip_labels, NA)

# Plot the tree without labels
plot(tree, main = "Phylogenetic Tree", cex = 0.4, show.tip.label = FALSE)

# Add colored labels for selected ASVs
p <- for (i in seq_along(tip_labels)) {
  if (tip_labels[i] %in% highlighted_ASVs_red) {
    tiplabels(text = tip_labels[i], tip = i, col = "red", font = 2, cex = 1.2)
  } else if (tip_labels[i] %in% highlighted_ASVs_blue) {
    tiplabels(text = tip_labels[i], tip = i, col = "blue", font = 2, cex = 1.2)
  }
}


#Add sequence simmularity colouring
# Install the seqinr package if you haven't done so
install.packages("seqinr")

# Load the seqinr package
library(seqinr)

# Specify the path for the FASTA file
fasta_file_path <- "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/flu_aasequence.fasta"

# Write sequences to the FASTA file
write.fasta(sequences = as.list(merged_data$Sequence), 
            names = merged_data$OTU, 
            file.out = fasta_file_path)

#Rhodo fasta
write.fasta(sequences = as.list(rhodobacteraceae_data$Sequence), 
            names = rhodobacteraceae_data$OTU, 
            file.out = fasta_file_path)

# Create the initial plot with the tree
p <- ggtree(tree) + 
  geom_tiplab(size = 2)

# Loop through each unique Order to shade clusters
# unique_orders <- unique(merged_data$Order)
# 
# 
# for (order in unique_orders) {
#   # Get the ASVs for the current Order
#   asvs_in_order <- merged_data$OTU[merged_data$Order == order]
#   
#   # Check if the ASVs are in the tree
#   asvs_in_order <- asvs_in_order[asvs_in_order %in% tree$tip.label]
#   
#   if (length(asvs_in_order) > 1) {  # Only consider orders with more than one ASV
#     # Identify the tip indices for the ASVs in this order
#     tip_indices <- which(tree$tip.label %in% asvs_in_order)
#     
#     # Find the MRCA node for these ASVs
#     mrca_node <- getMRCA(tree, tip_indices)
#     
#     # Highlight the MRCA node with a unique color
#     if (!is.null(mrca_node)) {
#       p <- p + geom_hilight(node = mrca_node, fill = sample(colors(), 1), alpha = 0.3)
#     }
#   }
# }

# Create a data frame for tip labels with colors based on selection
tip_labels <- data.frame(
  label = tree$tip.label,
  color = ifelse(tree$tip.label %in% highlighted_ASVs_red, "red",
                 ifelse(tree$tip.label %in% highlighted_ASVs_blue, "blue", "black"))
)

# Initialize the ggtree plot with custom tip labels and colors
p <- ggtree(tree) + 
  geom_tiplab(aes(label = label, color = color), 
              data = tip_labels, 
              size = 2) +  # Adjust label size as needed
  scale_color_identity()  # Use specified colors directly

# Now, plot MSA on top of the tree using the FASTA file path
msa_plot <- msaplot(p = p, fasta = fasta_file_path)

# Display the combined plot with a title
print(msa_plot + ggtitle("Phylogenetic Tree with MSA Overlay"))

# Plot MSA on top of the tree
# Now you can call msaplot() with the FASTA file path
msa_plot <- msaplot(p = p, fasta = fasta_file_path)


# Display the combined plot
print(msa_plot + ggtitle("Phylogenetic Tree with MSA Overlay"))

# Initialize the ggtree plot
p <- ggtree(tree) + 
  geom_tiplab(aes(label = tree$tip.label, 
                  color = ifelse(tree$tip.label %in% highlighted_ASVs_red, "red",
                                 ifelse(tree$tip.label %in% highlighted_ASVs_blue, "blue", "black"))),
              size = 2) +  # Adjust label size as needed
  scale_color_identity()  # Use specified colors directly

# Plot MSA on top of the tree using the FASTA file path
msa_plot <- msaplot(p = p, fasta = fasta_file_path)

# Display the combined plot with a title
print(msa_plot + ggtitle("Phylogenetic Tree with MSA Overlay"))


#Roseobacter clade ----
# Assuming `dist_matrix` is your distance matrix and `phylo_tree` is your tree
#Defining as Roseo if >88% similarity 
#citation for 89% similarity: 
#Buchan A, Gonzalez JM, Moran MA. 2005. Overview of the marine Roseobacter lineage. Appl. Environ. Microbiol. 71:5665–5677. 10.1128/AEM.71.10.5665-5677.2005. 

# Step 1: Extract ASV558 (MU42022) or ASV230 (PB2023) sequence from the data frame

asv558_sequence <- rhodobacteraceae_data$Sequence[rhodobacteraceae_data$OTU == "ASV558"][1]

#asv230_sequence <- rhodobacteraceae_data$Sequence[rhodobacteraceae_data$OTU == "ASV230"][1]

# Install and load Biostrings package if not already installed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("Biostrings")
}
library(Biostrings)

rhodobacteraceae_data <- subset(merged_data, Family == "Rhodobacteraceae")
#Extract unique sequences for the Rhodobacteraceae ASVs
sequences <- unique(rhodobacteraceae_data$Sequence)
names(sequences) <- unique(rhodobacteraceae_data$OTU)  # Assign ASV names to the sequences

dna_sequences <- DNAStringSet(sequences)

#ref_seq <- dna_sequences[["ASV230"]]

ref_seq <- dna_sequences[["ASV558"]]

# Initialize a vector to store results
roseobacter_matches <- list()
for (i in names(dna_sequences)) {
  # Skip if the current sequence is ASV558 itself
  if (i == "ASV230") next
  
  # Perform global pairwise alignment between ASV558 and the current sequence
  alignment <- pairwiseAlignment(ref_seq, dna_sequences[[i]], type = "global")
  
  # Calculate percentage identity
  identity <- mean(pid(alignment), na.rm = TRUE)  # Remove NAs
  
  # If identity is not NA and similarity is 89% or higher, add to the list of Roseobacter matches
  if (!is.na(identity) && identity >= 89) {
    roseobacter_matches[[i]] <- dna_sequences[[i]]
  }
}

# Check the results
print(roseobacter_matches)

# View the names of sequences classified as Roseobacter
names(roseobacter_matches)

# Add ASV558 to the results if you want to include it explicitly
roseobacter_matches[["ASV558"]] <- ref_seq
#roseobacter_matches[["ASV230"]] <- ref_seq

# Check the results
print(roseobacter_matches)

# View the names of sequences classified as Roseobacter
names(roseobacter_matches)

# Step 1: Filter the dataframe for ASVs in roseobacter_matches
roseobacter_df <- rhodobacteraceae_data[rhodobacteraceae_data$OTU %in% names(roseobacter_matches), ]

# Step 2: Extract the DNA sequences of the filtered ASVs
# Ensure that we’re extracting unique ASV sequences
roseobacter_sequences <- unique(roseobacter_df[c("OTU", "Sequence")])

sequences <- unique(roseobacter_sequences$Sequence)
names(sequences) <- unique(roseobacter_sequences$OTU)

dna_sequences <- DNAStringSet(sequences)

# Align the sequences
aligned_sequences <- AlignSeqs(dna_sequences) 
# Convert DNAStringSet to a character matrix
aligned_matrix <- as.matrix(aligned_sequences)

# Convert the matrix to a phyDat object
aligned_phyDat <- phyDat(aligned_matrix, type = "DNA")

# Check the phyDat object
print(aligned_phyDat)

# Compute the distance matrix
dist_matrix <- dist.ml(aligned_phyDat)

# Build the phylogenetic tree
tree <- NJ(dist_matrix)
plot(tree, main = "Roseobacter - Phylogenetic Tree", cex = 1)  # Adjust cex as needed

# highlighted_ASVs_red <- c("ASV11", "ASV88", "ASV198", "ASV178", "ASV201", "ASV471", "ASV613")
# highlighted_ASVs_blue <- c("ASV7", "ASV18")

# Get the tip labels
tip_labels <- tree$tip.label

# Add red and blue bold labels for the specified ASVs
for (i in seq_along(tip_labels)) {
  if (tip_labels[i] %in% highlighted_ASVs_red) {
    tiplabels(text = tip_labels[i], tip = i, col = "red", font = 2, cex = 1.5)
  } else if (tip_labels[i] %in% highlighted_ASVs_blue) {
    tiplabels(text = tip_labels[i], tip = i, col = "blue", font = 2, cex = 1.5)
  }
}


#Cirucular tree
library(ape)      # For plotting phylogenetic trees
library(phytools) # For additional phylogenetic functions

# Assuming 'tree' is your phylogenetic tree object
highlighted_ASVs_red <- c("ASV11", "ASV88", "ASV198", "ASV178", "ASV201", "ASV471", "ASV613")
highlighted_ASVs_blue <- c("ASV7", "ASV18")

# Plot the tree in a circular format
plot(tree, main = "Roseobacter - Circular Phylogenetic Tree", cex = 1, type = "fan")  # Use type = "fan"

# Get the tip labels
tip_labels <- tree$tip.label

# Add labels with colors for highlighted ASVs
for (i in seq_along(tip_labels)) {
  if (tip_labels[i] %in% highlighted_ASVs_red) {
    tiplabels(text = tip_labels[i], tip = i, col = "red", font = 2, cex = 1.3)
  } else if (tip_labels[i] %in% highlighted_ASVs_blue) {
    tiplabels(text = tip_labels[i], tip = i, col = "blue", font = 2, cex = 1.3)
  }
}

#add genus shading
genus_colors <- setNames(rhodobacteraceae_data$OTU, rhodobacteraceae_data$Genus)

# Plot the circular tree
plot(tree, main = "Roseobacter - Circular Phylogenetic Tree", cex = 0.5)

# Get the tip labels
tip_labels <- tree$tip.label

# Add labels with colors for highlighted ASVs
for (i in seq_along(tip_labels)) {
  if (tip_labels[i] %in% highlighted_ASVs_red) {
    tiplabels(text = tip_labels[i], tip = i, col = "red", font = 2, cex = 0.8)
  } else if (tip_labels[i] %in% highlighted_ASVs_blue) {
    tiplabels(text = tip_labels[i], tip = i, col = "blue", font = 2, cex = 0.8)
  } else {
    # Extract genus from the tip label
    genus <- rhodobacteraceae_data$Genus[rhodobacteraceae_data$ASV == tip_labels[i]]
    
    # Determine color based on genus
    color <- ifelse(genus %in% names(genus_colors), genus_colors[genus], "black")
    
    tiplabels(text = tip_labels[i], tip = i, col = color, font = 2, cex = 0.8)
  }
}

#Heatmaps ----
#Create a heatmap of Roseobacter abundance
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)

#create roseo dataframe
pseq <- PB2023_spat_not_rarefied_normalized
pseq_filt <- subset_samples(pseq, !Treatment %in% c("Continuous Probiotics", "James"))
any(taxa_sums(pseq_filt) == 0)
# pseq_filtered <- prune_taxa(taxa_sums(pseq_filt) > 0, pseq_filt)
# any(taxa_sums(pseq_filtered) == 0)
# pseq_filt <- pseq_filtered
pseq <- microbiome::transform(pseq_filt, "compositional")
pseq <- psmelt(pseq)
merged_data <- merge(pseq, sequence_ASVname_mapping_SMK, by.x = "OTU", by.y = "ASV", all.x = TRUE)
roseobacter_df <- merged_data[merged_data$OTU %in% names(roseobacter_matches), ]



# Filter and aggregate data for Roseobacter ASVs at the spat stage
roseobacter_spat_summary <- roseobacter_df %>%
  group_by(OTU, Treatment, Abundance) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = 'drop')


# Assuming roseobacter_spat_summary is your current data frame
heatmap_data <- roseobacter_spat_summary %>%
  # Summing duplicate values for OTU and Treatment pairs
  group_by(OTU, Treatment) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  # Reshaping the data to have OTU as rows and Treatment as columns
  pivot_wider(names_from = Treatment, values_from = Abundance, values_fill = 0)

# Convert to matrix for heatmap
roseobacter_matrix <- as.matrix(heatmap_data[, -1])
rownames(roseobacter_matrix) <- heatmap_data$OTU 


library(RColorBrewer)
pheatmap(
  roseobacter_matrix,
  color = colorRampPalette(c("white", "red", "darkred"))(50),
  border_color = "black",
  main = "",
  cluster_rows = TRUE, 
  cluster_cols = FALSE,
  show_rownames = TRUE,
  angle_col = 0,
  fontsize_col = 18,
  fontsize_row = 13,
  clustering_distance_rows = "manhattan"
)

# Perform hierarchical clustering
row_hclust <- hclust(dist(roseobacter_matrix), method = "complete")

# Cut the tree into clusters (e.g., 5 clusters)
k <- 3  # Number of clusters
row_clusters <- cutree(row_hclust, k = k)

# Create a data frame for row annotations
row_annotation <- data.frame(Cluster = factor(row_clusters))
rownames(row_annotation) <- rownames(roseobacter_matrix)

# Assign colors for the clusters
annotation_colors <- list(Cluster = c(
  "1" = "yellow",
  "2" = "orange",
  "3" = "forestgreen",
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
  fontsize_col = 13
)

#Boxplots ----

# Calculate average and standard deviation of Abundance for each OTU across treatment groups
roseobacter_stats <- roseobacter_df %>%
  group_by(OTU, Treatment) %>%
  summarise(
    Average_Abundance = mean(Abundance, na.rm = TRUE),
    Std_Abundance = sd(Abundance, na.rm = TRUE),
    .groups = 'drop'  # Drop grouping after summarising
  )


library(RColorBrewer)

library(viridis)  # For gradient color palettes

# Create a stacked bar plot of ASV abundance by Treatment with no gridlines and a gradient color palette
ggplot(roseobacter_stats, aes(x = Treatment, y = Average_Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack") +  # Change position to "stack"
  labs(title = "",
       x = "",  
       y = "Average Relative Abundance") +
  scale_fill_viridis_d(option = "F", direction = 1) +  # Use Viridis palette; adjust option for different colors
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=16, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#scatter plot
ggplot(roseobacter_stats, aes(x = Treatment, y = Average_Abundance, color = OTU)) +
  geom_jitter(width = 0.2, height = 0, size = 7, alpha = 0.7) + 
  geom_text(aes(label = OTU), vjust = -0.5, size = 4) + 
  labs(title = "",
       x = "",  
       y = "Average Relative Abundance") +
  scale_color_viridis_d(option = "F", direction = 1) +  # Use Viridis palette for colors
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=16, face = "bold"),
        panel.grid.major = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank())



#Testing roseobacter ----

pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq3 <- microbiome::transform(pseq, "compositional")

pseq3 <- psmelt(pseq3) 

roseobacter_df <- pseq3[pseq3$OTU %in% names(roseobacter_matches), ]

head(roseobacter_df)

# One-way ANOVA
anova_result <- aov(Abundance ~ Treatment, data = roseobacter_df)
summary(anova_result)
plot(anova_result)

#Many zeros
library(MASS)  # For negative binomial GLM
library(lme4)  # For mixed models

# GLM with negative binomial distribution
nb_glm <- glm.nb(Abundance ~ Treatment, data = roseobacter_df)
summary(nb_glm)

# GLMM with negative binomial distribution
nb_glmm <- glmer.nb(Abundance ~ Treatment + (1|Genetics), data = roseobacter_df)
summary(nb_glmm)

# Simulate residuals
qqnorm(residuals(nb_glmm))
plot(fitted(nb_glmm)~residuals(nb_glmm))
residuals_dharma <- simulateResiduals(fittedModel = nb_glmm)
plot(residuals_dharma)

#Try tweedie dist'd
library(mgcv)
tweedie_model <- mgcv::gam(Abundance ~ Treatment, 
                     family = tw(), 
                     data = roseobacter_df)
summary(tweedie_model)

qqnorm(residuals(tweedie_model))
residuals_dharma <- simulateResiduals(fittedModel = tweedie_model)
plot(residuals_dharma)

#manyglm

pseq2 <- MU42022_filtered_Oct92024

pseq2 <- subset_samples(pseq2, Age %in% c("Spat"))
pseq2 <- subset_samples(pseq2, !Genetics %in% c("4"))

fact1 = sample_data(pseq2)
fact = as.matrix.data.frame(fact1)
fact = as.data.frame(fact)

roseobacter_taxa <- names(roseobacter_matches)
pseq2_subset <- prune_taxa(roseobacter_taxa, pseq2)
ASV_data_cleaned <- as.data.frame(pseq2_subset@otu_table)
#rownames(ASV_data_cleaned) <- NULL


#many zeros - try method
# install.packages("zCompositions")
# library(zCompositions)
#roseobacter_otu_matrix[is.na(roseobacter_otu_matrix)] <- 0
#roseobacter_otu_czm <- cmultRepl(roseobacter_otu_matrix[, -1], method = "CZM", label = 0)
# too many zeros - does not work, could try on full dataset
#roseobacter_otu_bayesian <- cmultRepl(roseobacter_otu_matrix[, -1], method = "GBM", label = 0)
# roseobacter_otu_matrix[, -1] <- apply(roseobacter_otu_matrix[, -1], 2, function(x) as.numeric(as.character(x)))
# roseobacter_otu_bayesian <- cmultRepl(roseobacter_otu_matrix[, -1], method = "GBM", label = 0)

dat_mvabund <- mvabund(ASV_data_cleaned)

library(mvabund)

Mod = manyglm(dat_mvabund ~ Treatment * Genetics, family="negative.binomial", data = fact, composition=TRUE, show.warning = TRUE)
#Warning: singular matrix in betaEst: An eps*I is added to the singular matrix.
plot(Mod)
summary(Mod)

#Using CZM method on pseq object
library(compositions)

# Extract OTU table and convert to a matrix
otu_matrix <- otu_table(pseq2)

# Apply the CZM method
otu_czm <- cmultRepl(as.matrix(otu_matrix), method = "CZM", label = 0)

#centered log-normalization
otu_czm_clr <- clr(otu_czm)

# Convert clr data to a matrix
otu_czm_clr_matrix <- as.matrix(otu_czm_clr)

# Ensure OTU names are consistent
rownames(otu_czm_clr_matrix) <- rownames(otu_table(pseq2))

# Create the new OTU table
otu_table_new <- otu_table(otu_czm_clr_matrix, taxa_are_cols = TRUE)

# Recreate the phyloseq object
physeq_new <- phyloseq(otu_table_new, sample_data(pseq2), tax_table(pseq2))

# Save the data frame as a tab-delimited text file
write.table(roseobacter_df, file = "roseobacter_df_PB2023.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Pivot the data frame to wide format
roseobacter_wide_df <- roseobacter_df %>%
  pivot_wider(names_from = OTU, values_from = Abundance, values_fill = list(Abundance = 0))

# View the resulting wide-format data frame
View(roseobacter_wide_df)


#test if avergae roseobacter abundance is different between treatments
library(dplyr)

# Summarize average abundance by Sample and Treatment
sample_abundance <- roseobacter_df %>%
  filter(!Treatment %in% c("Continuous Probiotics", "James")) %>%
  group_by(Sample.ID, Treatment) %>%
  summarize(avg_abundance = mean(Abundance, na.rm = TRUE))

library(lme4)



