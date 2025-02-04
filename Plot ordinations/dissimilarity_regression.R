#Dissimilaritiy comparisons

#Phylogenetic relatedness between ASVs
library("devtools")
library(phyloseq)
library(microbiome)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(ggtree)
library(ape)
library(vegan)
library(reshape2)
library(dplyr)

MU42022_filtered_Oct92024 <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/MU42022_filtered_Oct92024.rds")

pseq <- MU42022_filtered_Oct92024

pseq <- microbiome::transform(pseq, "compositional")

#remove PB ASVs
otu_table_data <- otu_table(pseq)

# Subset the OTU table to exclude ASV7 and ASV18 columns
otu_table_filtered <- otu_table_data[, !(colnames(otu_table_data) %in% c("ASV7", "ASV18"))]

# Reconstruct the phyloseq object with the filtered OTU table
pseq_filtered <- phyloseq(otu_table(otu_table_filtered, taxa_are_rows = FALSE), sample_data(pseq), tax_table(pseq))

# Verify the changes
pseq_filtered
pseq <- pseq_filtered


otu_matrix <- otu_table(pseq)
bray_dist <- vegdist(otu_matrix, method = "bray")

# Step 2: Convert the distance matrix to a data frame
bray_df <- as.matrix(bray_dist)
bray_df <- melt(bray_df)  # Converts to long format
colnames(bray_df) <- c("Sample1", "Sample2", "Dissimilarity")


# Step 3: Merge with metadata
metadata <- data.frame(sample_data(pseq))  # Extract sample metadata
bray_df <- bray_df %>%
  left_join(metadata, by = c("Sample1" = "Sample.ID")) %>%
  rename(Treatment1 = Treatment, Time1 = Age) %>%  # Adjust if necessary
  left_join(metadata, by = c("Sample2" = "Sample.ID")) %>%
  rename(Treatment2 = Treatment, Time2 = Age)

# Step 4: Filter or summarize as needed for comparison between Treatments and Times
bray_df_filtered <- bray_df %>%
  filter(Treatment1 != Treatment2)  # Compare only between different treatments
# You can also add conditions to compare specific time points, if needed.

# Step 5: Create dot plot
ggplot(bray_df_filtered, aes(x = Treatment1, y = Dissimilarity, color = Treatment2)) +
  geom_point(alpha = 0.6, size = 3) +
  facet_wrap(~ Time1, scales = "free_y") +
  theme_minimal() +
  labs(title = "Dissimilarities Between Treatments Over Time",
       x = "Treatment (Sample 1)", y = "Dissimilarity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



bray_df_filtered <- bray_df %>%
  filter(Treatment1 != Treatment2) %>%
  filter(Time1 == Time2)  # Keep only comparisons with the same Age (Time)

# Plot using ggplot
ggplot(bray_df_filtered, aes(x = Treatment1, y = Dissimilarity, color = Treatment2)) +
  geom_point(aes(shape = factor(Time1)), alpha = 0.7, size = 3) +  # Shape represents Age
  facet_wrap(~ Time1, scales = "free_y") +  # Create separate plots for each time point (Age)
  theme_minimal() +
  labs(title = "Dissimilarity Between Treatments at Same Age",
       x = "Treatment 1", y = "Dissimilarity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Filter for comparisons between Control and Probiotics treatments at the same Age
bray_df_filtered <- bray_df %>%
  filter((Treatment1 == "Probiotics" & Treatment2 == "Probiotics + HT") | 
           (Treatment1 == "Probiotics + HT" & Treatment2 == "Probiotics")) %>%
  filter(Time1 == Time2)  # Keep only comparisons with the same Age (Time)

bray_df_filtered <- bray_df %>%
  filter((Treatment1 == "Probiotics" & Treatment2 == "Control") | 
           (Treatment1 == "Control" & Treatment2 == "Probiotics")) %>%
  filter(Time1 == Time2)  # Keep only comparisons with the same Age (Time)

# Plot using ggplot
ggplot(bray_df_filtered, aes(x = Time1, y = Dissimilarity, color = Treatment2)) +
  geom_point(alpha = 0.7, size = 3) +  # Shape represents Age
  #facet_wrap(~ Time1, scales = "free_y") +  # Create separate plots for each time point (Age)
  theme_minimal() +
  labs(title = "Dissimilarity Between Control and Probiotics Treatments at Same Age",
       x = "Time (Age)", y = "Dissimilarity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels


#one way comparison between treatments ----
bray_df_filtered <- bray_df %>%
  filter((Treatment1 == "Probiotics" & Treatment2 == "Control")) %>%
  filter(Time1 == Time2)  # Keep only comparisons with the same Age (Time)

#one way comparison
bray_df_filtered <- bray_df %>%
  filter((Treatment1 == "Probiotics" & Treatment2 == "Control")) %>%
  filter(Time1 == Time2)  # Keep only comparisons with the same Age (Time)


bray_df_filtered <- bray_df %>%
  filter((Treatment1 == "Probiotics + HT" & Treatment2 == "Control")) %>%
  filter(Time1 == Time2)  # Keep only comparisons with the same Age (Time)

bray_df_filtered <- bray_df %>%
  filter((Treatment1 == "Probiotics + HT" & Treatment2 == "Probiotics")) %>%
  filter(Time1 == Time2)  # Keep only comparisons with the same Age (Time)



lm_model <- lm(Dissimilarity ~ Time1, data = bray_df_filtered)
summary(lm_model)
plot(lm_model)

# Extract coefficients and R²
equation <- paste("y = ", round(coef(lm_model)[2], 2), "*x + ", round(coef(lm_model)[1], 2))
r_squared <- summary(lm_model)$r.squared
r_squared_text <- paste("R² = ", round(r_squared, 2))

library(ggpmisc)

# Create the plot with regression line and equation
ggplot(bray_df_filtered, aes(x = Time1, y = Dissimilarity)) +
  geom_point(alpha = 0.7, size = 3, color = "orange") +  # Plot points with limegreen color
  geom_smooth(method = "lm", se = FALSE, aes(group = Treatment2), color = "orange") +  # Add regression line with limegreen color
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(adj.r.squared), sep = "*\", \"*")),
               formula = y ~ x, parse = TRUE) +  # Display equation on the plot
  annotate("text", x = 1.5, y = 0, label = paste(equation, r_squared_text), size = 4, color = "black") +  # Add equation and R²
  theme_bw() +
  labs(title = "Dissimilarity Between C/PBH",
       x = "", y = "Dissimilarity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()) +  # Remove minor gridlines
  ylim(0, 0.8)


