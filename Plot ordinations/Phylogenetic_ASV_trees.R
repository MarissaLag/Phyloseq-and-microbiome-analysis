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

pseq2 <- MU42022_filtered_Oct92024
pseq2 <- psmelt(pseq2) 
head(pseq2)

#Merge ASV seq with metadata
#'MU42022_sequence_ASVname_mapping'
colnames(MU42022_sequence_ASVname_mapping) <- c("ASV", "Sequence")

head(pseq2@tax_table)
head(MU42022_sequence_ASVname_mapping)

# Merge the ASV table with the sequence data based on ASV number
merged_data <- merge(pseq2, MU42022_sequence_ASVname_mapping, by.x = "OTU", by.y = "ASV", all.x = TRUE)
View(merged_data)
str(merged_data)

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
plot(tree, main = "Phylogenetic Tree", cex = 0.35, show.tip.label = FALSE)

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

# Step 1: Extract ASV558 sequence from the data frame
asv558_sequence <- rhodobacteraceae_data$Sequence[rhodobacteraceae_data$OTU == "ASV558"][1]

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

ref_seq <- dna_sequences[["ASV558"]]

# Initialize a vector to store results
roseobacter_matches <- list()

for (i in names(dna_sequences)) {
  # Skip if the current sequence is ASV558 itself
  if (i == "ASV558") next
  
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
plot(tree, main = "Roseobacter - Phylogenetic Tree", cex = 0.5)  # Adjust cex as needed

highlighted_ASVs_red <- c("ASV11", "ASV88", "ASV198", "ASV178", "ASV201", "ASV471", "ASV613")
highlighted_ASVs_blue <- c("ASV7", "ASV18")

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


#Cirucular tree
library(ape)      # For plotting phylogenetic trees
library(phytools) # For additional phylogenetic functions

# Assuming 'tree' is your phylogenetic tree object
highlighted_ASVs_red <- c("ASV11", "ASV88", "ASV198", "ASV178", "ASV201", "ASV471", "ASV613")
highlighted_ASVs_blue <- c("ASV7", "ASV18")

# Plot the tree in a circular format
plot(tree, main = "Roseobacter - Circular Phylogenetic Tree", cex = 0.5)  # Use type = "fan"

# Get the tip labels
tip_labels <- tree$tip.label

# Add labels with colors for highlighted ASVs
for (i in seq_along(tip_labels)) {
  if (tip_labels[i] %in% highlighted_ASVs_red) {
    tiplabels(text = tip_labels[i], tip = i, col = "red", font = 2, cex = 0.8)
  } else if (tip_labels[i] %in% highlighted_ASVs_blue) {
    tiplabels(text = tip_labels[i], tip = i, col = "blue", font = 2, cex = 0.8)
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

# Filter and aggregate data for Roseobacter ASVs at the spat stage
roseobacter_spat_summary <- roseobacter_df %>%
  filter(OTU %in% names(roseobacter_matches), Age == "Spat", !Genetics == "4") %>%
  select(OTU, Treatment, Abundance) %>%
  group_by(OTU, Treatment) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = 'drop') %>%
  spread(key = Treatment, value = Total_Abundance, fill = 0)

#Remove grouping and calculate total abundance across treatments
roseobacter_spat_summary <- roseobacter_spat_summary %>%
  ungroup() %>%  # Remove grouping to ensure numeric-only columns for rowSums
  mutate(Total = rowSums(select(., -OTU))) %>%  # Calculate total abundance across treatments
  filter(Total > 0) %>%  # Filter out ASVs with zero abundance across treatments
  select(-Total)  # Remove the Total column after filtering

# Convert to matrix for heatmap
roseobacter_matrix <- as.matrix(roseobacter_spat_summary[,-1])
rownames(roseobacter_matrix) <- roseobacter_spat_summary$OTU

library(RColorBrewer)
pheatmap(
  roseobacter_matrix,
  border_color = "black",
  color = brewer.pal(9, "Reds"),
  main = "",
  cluster_rows = TRUE, 
  cluster_cols = FALSE,
  show_rownames = FALSE,
  angle_col = 0,
  fontsize_col = 12,
  clustering_distance_rows = "manhattan"
)

#relative abundance heatmap
# Summarize Roseobacter ASV abundance by Treatment at the Spat stage
roseobacter_spat_summary <- roseobacter_df %>%
  filter(OTU %in% names(roseobacter_matches), Age == "Spat", !Genetics == "4") %>%
  select(OTU, Treatment, Abundance) %>%
  group_by(OTU, Treatment) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = 'drop') %>%
  spread(key = Treatment, value = Total_Abundance, fill = 0)  # Fill missing values with 0

# Filter out ASVs with zero total abundance across all treatments
roseobacter_spat_summary <- roseobacter_spat_summary %>%
  rowwise() %>%
  filter(sum(c_across(-OTU)) > 0) %>%  # Keep only rows with non-zero total abundance
  ungroup()

# Calculate relative abundance
roseobacter_spat_summary <- roseobacter_spat_summary %>%
  rowwise() %>%
  mutate(across(-OTU, ~ . / sum(c_across(-OTU)))) %>%  # Calculate relative abundance
  ungroup()

# Convert to matrix for heatmap
roseobacter_matrix <- as.matrix(roseobacter_spat_summary[,-1])
rownames(roseobacter_matrix) <- roseobacter_spat_summary$OTU

# Plot the heatmap with relative abundance
pheatmap(
  roseobacter_matrix,
  border_color = "black",
  color = brewer.pal(9, "Reds"),
  main = "",
  cluster_rows = TRUE, 
  cluster_cols = FALSE,
  show_rownames = FALSE,
  angle_col = 0,
  fontsize_col = 12,
  clustering_distance_rows = "manhattan"
)


#Boxplots ----
roseobacter_spat<- rhodobacteraceae_data %>%
  filter(OTU %in% names(roseobacter_matches), Age == "Spat", !Genetics == "4")

library(RColorBrewer)

library(viridis)  # For gradient color palettes

# Create a stacked bar plot of ASV abundance by Treatment with no gridlines and a gradient color palette
ggplot(roseobacter_spat, aes(x = Treatment, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack") +  # Change position to "stack"
  labs(title = "",
       x = "",  
       y = "Total Abundance") +
  scale_fill_viridis_d(option = "F", direction = -1) +  # Use Viridis palette; adjust option for different colors
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())




