###doing correlations and heatmaps from mb2021 data
install.packages("Hmisc")
library("devtools")
library(phyloseq)
library(microbiome)
library(Hmisc)

#load data

Marissa_Osyter <- readRDS("~/mb2021/Marissa_Osyter.rds")

pseq <- Marissa_mb2021_filtered_20240203

#create objects

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

# Calculate Spearman's correlation matrix
cor_matrix <- rcorr(as.matrix(OTU), type = "spearman")

cor_matrix <- cor_matrix$r

heatmap(cor_matrix, 
        symm = TRUE,      # To show both sides of the diagonal (symmetrical)
        col = viridis::viridis(100),  # Color palette (you can choose any other)
        main = "Spearman's Correlation Heatmap",
        xlab = "Sample", ylab = "Sample",
        cexCol = 0.7, cexRow = 0.7, # Adjust label size
        margins = c(10, 10))  # Adjust margins if needed

##pearson correlation heatmap

install.packages("ggplot2")
install.packages("reshape2")
install.packages("viridis")

library(ggplot2)
library(reshape2)
library(viridis)

# Convert the correlation matrix to a tidy format
cor_data <- melt(cor_matrix)

# Create the correlation plot with Viridis color palette
ggplot(data = cor_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pearson's Correlation Heatmap",
       x = "Taxa", y = "Taxa", fill = "Correlation")


##run pearson's correlation and heatmap with only larval samples

(ggplot2)
library(reshape2)
library(viridis)

excluded_samples <- c("T10r1", "T10r2", "T10r3", "T11r1", "T11r3", "T12r1", "T12r2", "T12r3", "T13r1","T13r2","T13r3", "T14r1", "T14r2", "T15r1", "T15r2", "T16r1", "T16r2", "T16r3", "T1r1", "T1r2", "T1r3", "T2r1", "T2r3", "T3r1", "T3r2", "T3r3", "T4r1", "T4r2", "T9r1", "T9r3", "T11r2", "T14r3", "T2r2")

microbiome_data_filtered <- OTU[, !colnames(OTU) %in% excluded_samples]

cor_matrix <- cor(microbiome_data_filtered, method = "pearson")

# Convert the correlation matrix to a tidy format
cor_data <- reshape2::melt(cor_matrix)

# Plot the correlation heatmap
ggplot(data = cor_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pearson's Correlation Heatmap",
       x = "Taxa", y = "Taxa", fill = "Correlation")

# Create the correlation plot with custom axis labels
ggplot(data = cor_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pearson's Correlation Heatmap",
       x = "", y = "", fill = "Correlation") +  # Clear existing axis labels
  scale_x_discrete(labels = "Control, 1 DPF", "Control, 18 DPF", "Control, 3 DPF", "Control, 3 DPF", "HS, 1 DPF", "HS, 18 DPF", "HS, 3 DPF", "LS, 18 DPF", "Control, 1 DPF", ") +  # Set custom x-axis labels
  scale_y_discrete(labels = )    # Set custom y-axis labels
