#### Inval Tutorial Draft ----

### Source:  https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html 

#Packages ----

install.packages("indicspecies")
library(indicspecies)
#Set data ----

Marissa_MU42022_rare <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rare.rds")

pseq <- Marissa_MU42022_rarefied_20231016
pseq <-`Filtered_Rarified_MU42022_23-12-13`

#Load objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree


#Extract abundance matrix ----
#from the phyloseq object using phyloseq

OTU1 = as(OTU, "matrix")
Metadata = as(Metadata, "matrix")
write.csv(OTU1, file="Data_fram_1_James.cvs",row.names=TRUE)
write.csv(Metadata, file="Data_fram_1_James_Meta.cvs",row.names=TRUE )

write.table(OTU1,file="data_table_James.csv",sep=",",dec = " ")
write.table(Metadata,file="data_table_James_Meta.csv",sep=",",dec = " ")
####Format to example data and reload below for actual test 

#reload edited table

pc_FUN = read.csv("data_table_James.csv", header= TRUE)

View(pc_FUN)

####Test ASVs ----

#Inverse data
funi_df<- t(pc_FUN)

View(funi_df)
View(pc_FUN)

#if need to remove samples

rows_to_remove <- 1
new_matrix <- pc_FUN[-rows_to_remove, ]
View(new_matrix)
pc_FUN <- new_matrix

#to subset James' samples into antibiotics and control only 

AntiB <- c("Antibiotics")

AntiB <- pc_FUN[pc_FUN[, 3] %in% AntiB, ]

View(AntiB)

Control <- c("Control")

Control <- pc_FUN[pc_FUN[, 3] %in% Control, ]

View(Control)

###make into a matrix and populate::: This tells r what is metadata and what is the actual data ... Below 5-952 are the coloumns that are the data

matrix_F = pc_FUN[ ,9:350]

matrix_F = AntiB[ ,9:350]

matrix_F = Control[ ,9:350]

# Select covariates to test
time_a_F <- pc_FUN$Treatment_Age

time_a_F = pc_FUN$Age

time_a_F = pc_FUN$Treatment

time_a_F = Control$Age

time_a_F = AntiB$Age

### Run test 
inv_F = multipatt(matrix_F, time_a_F, func = "r.g", control = how(nperm=999))
summary(inv_F)

#heat maps ----
library(pheatmap)
library(ggplot2)


subset_data <- inv_F$str[, selected_columns]
comb_matrix <- inv_F$comb
pheatmap(comb_matrix)

View(inv_F)

str(inv_F)

View(inv_F$str)

View(inv_F$sign)

#for some reason some rows say "NA"
inv_F$str <- inv_F$str[rownames(inv_F$str) != "ASV384", ]
cleaned_matrix <- inv_F$str[complete.cases(inv_F$str), ]
inv_F <-cleaned_matrix

#To view all ASVs and samples
pheatmap(inv_F$str)


# Assuming p.value is the column name for p-values in the sign section
significant_asvs_indices <- which(inv_F$sign$p.value < 0.01)  # Substitute the threshold as per your significance criterion
significant_asvs_matrix <- inv_F$str[significant_asvs_indices, ]
View(significant_asvs_matrix)
pheatmap(significant_asvs_matrix)

#exclude certain columns
selected_columns <- c("Day 01", "Day 03", "Day 06", "Day 15", "Spat")  # Replace with your actual column names
selected_columns <- c("Day 01", "Day 03", "Day 06", "Day 08", "Day 15")

# Create a subset of your data including only the selected columns
subset_data <- significant_asvs_matrix[, selected_columns]
subset_data <- inv_F$str[, selected_columns]

new_column_order <- c("Day 01", "Day 03", "Day 06", "Day 15", "Spat")
new_column_order <- c("Day 01", "Day 03", "Day 06", "Day 08", "Day 15")

# Create the heatmap using the subset of your data
pheatmap(subset_data, column_order = new_column_order, cluster_cols = FALSE, show_rownames = FALSE, main = "Control")
pheatmap(subset_data, column_order = new_column_order, cluster_cols = FALSE, show_rownames = FALSE, main = "Antibiotics")


#Creating annotations in a heatmap
#source: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/ 

install.packages("dendextend")
library(dendextend)

#create dendrogram of variable (e.g., genes or ASVs)

my_hclust_ASV <- hclust(dist(data_subset), method = "complete")

as.dendrogram(my_hclust_ASV) %>%
  plot(horiz = TRUE)

#cut tree into clusters - here splitting into 2 clusters

my_ASV_col <- cutree(tree = as.dendrogram(my_hclust_ASV), k = 2)

#rename to cluster name

my_ASV_col <- data.frame(cluster = ifelse(test = my_ASV_col == 1, yes = "cluster 1", no = "cluster 2"))

head(my_gene_col)

#add multiple row annotations to a heatmap and below I'll add some random annotations.
set.seed(1984)
my_random <- as.factor(sample(x = 1:2, size = nrow(my_ASV_col), replace = TRUE))
my_ASV_col$random <- my_random

head(my_ASV_col)

#add some column annotations and create the heatmap.

my_sample_col <- data.frame(sample = rep(c("Control", "Antibiotics"), c(4,2)))
row.names(my_sample_col) <- colnames(data_subset)
pheatmap(data_subset, annotation_row = my_gene_col, annotation_col = my_sample_col)