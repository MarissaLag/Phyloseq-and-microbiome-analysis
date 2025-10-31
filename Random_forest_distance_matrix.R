#Force directed network as done in paper "Microbiomes of the Sydney Rock Oyster are acquired through both vertical and horizontal transmission"

library(phyloseq)
library(dplyr)

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

pseq <- Sam_all_samples_partial_rare_CSS

# Transform to relative abundance
pseq.rel <- transform_sample_counts(pseq, function(x) x / sum(x))

# Melt data to long format
df <- psmelt(pseq.rel)

#Remove bad columns
df <- df[, !grepl("^X\\.[0-9]+$", names(df))]

edges <- df %>%
  select(OTU, Stage) %>%
  distinct()


library(igraph)

# Create a graph from the edge list
g <- graph_from_data_frame(edges, directed = FALSE)

# Detect communities (modules) in the graph
communities <- cluster_fast_greedy(g)

# Assign colors to modules
V(g)$color <- communities$membership

# Plot the graph
plot(g, vertex.size=5, vertex.label=NA, layout=layout_with_fr)

# Assign colors to modules
V(g)$color <- communities$membership

# Plot the graph
plot(g, vertex.size=5, vertex.label=NA, layout=layout_with_fr)



library(igraph)
library(dplyr)

# Step 1: Make sure your edge list connects ASVs to Sample Types
edges <- df %>%
  filter(Abundance > 0.001) %>%
  select(OTU, Stage) %>%
  distinct() %>%
  filter(!is.na(Stage))  # Remove NA types if present

# Step 2: Create graph
g <- graph_from_data_frame(edges, directed = FALSE)

# Step 3: Assign node types
V(g)$type <- ifelse(V(g)$name %in% edges$OTU, "ASV", "Sample.ID")

# Step 4: Assign colors
sample_types <- unique(edges$Stage)
sample_colors <- setNames(RColorBrewer::brewer.pal(length(sample_types), "Set3"), sample_types)

V(g)$color <- ifelse(
  V(g)$type == "Sample.ID",
  sample_colors[V(g)$name],
  "grey70"  # ASVs all grey, or use modularity if you prefer
)

# Optional: change shape
V(g)$shape <- ifelse(V(g)$type == "Sample.ID", "square", "circle")

# Step 5: Plot
plot(g,
     vertex.label = NA,
     vertex.size = 5,
     layout = layout_with_fr)

#Try assigning ASVs with certain sample types
# Step 1: Get ASV-Stage pairs
asv_stage <- df %>%
  select(OTU, Stage) %>%
  filter(!is.na(Stage)) %>%
  group_by(OTU, Stage) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(OTU) %>%
  slice_max(count, n = 1, with_ties = TRUE)

library(RColorBrewer)

stages <- unique(asv_stage$Stage)
stage_colors <- setNames(RColorBrewer::brewer.pal(length(stages), "Set3"), stages)





library(dplyr)
library(tidyr)
library(igraph)

# Filter to relevant ASVs (as per your earlier filtering)
df_net <- df %>% 
  filter(Abundance > 0) %>% 
  select(Sample.ID, OTU, Abundance, Stage)

# Create edge list: Sample nodes to ASV nodes
# edges with character names (Sample.ID and OTU)
edges_char <- df_net %>%
  filter(Abundance > 0) %>%
  select(Sample.ID, OTU) %>%
  distinct() # optional: unique edges only

# nodes data frame
samples <- data.frame(name = unique(edges_char$Sample.ID), type = "Sample")
asvs <- data.frame(name = unique(edges_char$OTU), type = "ASV")
nodes <- rbind(samples, asvs)

# Now create graph with character names in edges
g <- graph_from_data_frame(edges_char, vertices = nodes, directed = FALSE)

print(g)
summary(g)

# Assign colors for samples by their Stage
# Create a named vector of colors for stages present in your nodes
library(RColorBrewer)
sample_nodes <- V(g)[type == "Sample"]
sample_stages <- nodes$name[which(nodes$type == "Sample")] # vector of sample names
# You will need to join your sample metadata to get stages per sample name:
sample_metadata <- df_net %>% distinct(Sample.ID, Stage)
stage_colors <- setNames(brewer.pal(length(unique(sample_metadata$Stage)), "Set3"),
                         unique(sample_metadata$Stage))

# Create a vector matching each sample vertex name to its stage color
V(g)$color <- ifelse(V(g)$type == "Sample",
                     stage_colors[sample_metadata$Stage[match(V(g)$name, sample_metadata$Sample.ID)]],
                     "grey70") # ASVs colored grey

# Assign shapes by type
V(g)$shape <- ifelse(V(g)$type == "Sample", "square", "circle")

library(igraph)

comm <- cluster_louvain(g)
V(g)$community <- comm$membership

plot(g,
     vertex.label = NA,
     vertex.size = 5,
     layout = layout_with_fr,
     main = "Sample-ASV network with modularity communities")

install.packages("randomForest")
library(randomForest)

# Make ASV abundance matrix
feature_table <- df_net %>%
  select(Sample.ID, OTU, Abundance) %>%
  pivot_wider(names_from = OTU, values_from = Abundance, values_fill = 0)

# Make sure Stage is included per Sample.ID
sample_stages_rf <- df_net %>%
  distinct(Sample.ID, Stage)

rf_data <- feature_table %>%
  left_join(sample_stages_rf, by = "Sample.ID")

# Remove rows with NA in Stage
rf_data_clean <- rf_data %>% filter(!is.na(Stage))

# Extract features and labels again from the cleaned data
rf_features_clean <- rf_data_clean %>% select(-Sample.ID, -Stage)
rf_labels_clean <- factor(rf_data_clean$Stage)

# Now run randomForest
set.seed(123)
rf_model <- randomForest(x = rf_features_clean, y = rf_labels_clean, ntree = 10001)

print(rf_model)

#List important ASVs
importance_vals <- importance(rf_model)
varImpPlot(rf_model)  # Visualize overall variable importance

# Extract important matrix
importance_vals <- importance(rf_model)

# Sort ASVs by MeanDecreaseGini (or MeanDecreaseAccuracy)
top_asvs <- sort(importance_vals[, "MeanDecreaseGini"], decreasing = TRUE)

# Print top 10 important ASVs
head(top_asvs, 10)

# Plot variable importance
varImpPlot(rf_model, n.var = 20, main = "Top 20 Important ASVs")

install.packages("caret")
library(caret)

# Prepare data frame for caret
rf_data <- data.frame(rf_features_clean)
rf_data$Stage <- rf_labels_clean

# Define train control for LOOCV
train_control <- trainControl(method = "LOOCV")

# Train model using caret with LOOCV
set.seed(42)
rf_caret <- train(Stage ~ ., data = rf_data,
                  method = "rf",
                  trControl = train_control,
                  ntree = 1001)

# View results
print(rf_caret)

#Heat map of top ASVs
library(pheatmap)

# Subset samples by stage or order by stage
sample_order <- order(sample_metadata$Stage)
abundance_ordered <- abundance_matrix[, sample_order]

abundance_matrix <- df %>%
  filter(!is.na(Sample.ID), !is.na(OTU)) %>%
  select(OTU, Sample.ID, Abundance, Stage) %>%
  pivot_wider(names_from = Sample.ID, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("OTU") %>%
  as.matrix()

# Create annotation for heatmap columns
annotation_col <- data.frame(Stage = sample_metadata$Stage[sample_order])
rownames(annotation_col) <- colnames(abundance_ordered)

# Plot heatmap
pheatmap(abundance_ordered,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "ASV abundance heatmap by Stage")

library(igraph)

# Create graph
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# Assign colors for stages
stage_colors <- setNames(RColorBrewer::brewer.pal(length(unique(nodes$Stage)), "Set3"),
                         unique(nodes$Stage))

# Color vertices: samples by Stage color, ASVs grey or by importance
V(g)$color <- ifelse(V(g)$type == "Sample",
                     stage_colors[V(g)$Stage],
                     ifelse(V(g)$name %in% names(top_asvs[1:20]), "red", "grey70"))

# Set shapes
V(g)$shape <- ifelse(V(g)$type == "Sample", "square", "circle")

# Plot with force-directed layout
plot(g,
     vertex.label = NA,
     vertex.size = 5,
     layout = layout_with_fr)

