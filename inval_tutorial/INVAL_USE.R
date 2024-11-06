#### Inval Tutorial Draft ----

### Source:  https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html 


#to look for differentially abundant ASVs using indicator species analysis


#Packages ----

install.packages("indicspecies")
library(indicspecies)
library(ggplot2)
library(dplyr)
library(microbiome)
library(phyloseq)
library(RColorBrewer)

#Set data ----

Marissa_MU42022_rare <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rare.rds")

pseq <- Marissa_MU42022_rare_nochloro

pseq <- MU42022_filtered_NOT_rarefied #579 taxa

#pseq <- Marissa_mb2021_filtered_20240203

pseq <- mb2021_filtered_NOT_rarefied #1007 taxa

#Load objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree


#Extract abundance matrix ----
#from the phyloseq object using phyloseq

OTU1 = as(OTU, "matrix")
write.csv(OTU1, file="Data_fram_1.cvs",row.names=TRUE)

write.table(OTU1,file="data_table_PB2023_unrarefied.csv",sep=",",dec = " ")



####Format to example data and reload below for actual test 

#reload edited table
data_table <- read.csv("data_table_PB2023_unrarefied.csv")

pc_FUN = read.csv("data_table_PB2023_unrarefied.csv", header= TRUE)

#pc_FUN <- data_table_mb2021_unrarefied

pc_FUN <- data_table_MU_2022_rarefied

#if removing samples ----

#for MU42022, remove PB on day 01
pc_FUN <- pc_FUN[, !colnames(pc_FUN) %in% c("ASV7", "ASV18")]

#for mb2021 project remove remaining day 3 samples and T9 spat data

pc_FUN <- pc_FUN[!pc_FUN$`Time-point` == "3 dpf", ]


#for mb2021 project remove tank 9 from pc_Fun for mb2021
pc_FUN <- pc_FUN[!pc_FUN$`Treatment` == "James", ] 
pc_FUN <- pc_FUN[!pc_FUN$`Treatment` == "Continuous-Probiotics", ] 

View(pc_FUN)

#Day 1 only 

pc_FUN <- pc_FUN[pc_FUN$'Age' == "Day 01", ]

#Spat only

pc_FUN <- pc_FUN[pc_FUN$`Age` == "Spat", ]

#If present, filter NAs in column 1

pc_FUN <- pc_FUN[-1, ]

  
####Test ASVs ----

#Inverse data
funi_df<- t(pc_FUN) 

###make into a matrix and populate::: This tells r what is metadata and what is the actual data ... Below 5-952 are the coloumns that are the data

#Note: sum columns add up to zero so you may get an error
#matrix_F = pc_FUN[ ,6:1012] 

#mb2021
matrix_F = pc_FUN[ ,7:190]

### Make the equation. Saying we want to examine specific column of metadata
time_a_F = pc_FUN$Treatment

### Run test 
inv_F_spat = multipatt(matrix_F, time_a_F, func = "r.g", control = how(nperm=9999))
results <- summary(inv_F_spat)


#save results
write.csv(inv_F, "Spat_INVALsummary_results.csv", row.names = TRUE)


#Bubble plot results ----
#MU42022 project
#ASV7 and 18 unique to PB/PBH treatments
#plotting ASV 7 and 18 on day 1 data as was differentially abundant on day 1 (p < 0.05)


pc_FUN <- pc_FUN[pc_FUN$Age == "Day 01", ]


ggplot(pc_FUN, aes(x=Treatment, y= ASV, size = Abundance, color = Treatment)) + 
  geom_point(alpha=0.7)+ 
  scale_size(range = c(.1, 10)) +
  scale_colour_ipsum() +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("") +
  xlab("") +
  theme(legend.position = "none")

#plotting ASV 376 and 88 on spat data as was differentially abundant in spat (p < 0.05)

pc_FUN <- pc_FUN[pc_FUN$Age == "Spat", ]

pc_FUN <- pc_FUN[pc_FUN$ASV %in% c("ASV88", "ASV376"), ]


ggplot(pc_FUN, aes(x=Treatment, y= Genus, size = Abundance, color = Treatment)) + 
  geom_point(alpha=0.7)+ 
  scale_size(range = c(.1, 10)) +
  scale_colour_ipsum() +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("") +
  xlab("") +
  theme(legend.position = "none")

#heat maps ----
library(pheatmap)
library(ggplot2)


subset_data <- inv_F$str[, selected_columns]
comb_matrix <- inv_F$comb
pheatmap(comb_matrix)

View(inv_F)

str(inv_F)

View(inv_F$str)

View(inv_F_day1$sign)

#for some reason some rows say "NA"
inv_F$str <- inv_F$str[rownames(inv_F$str) != "ASV384", ]
cleaned_matrix <- inv_F$str[complete.cases(inv_F$str), ]
inv_F <-cleaned_matrix

#To view all ASVs and samples
pheatmap(inv_F_spat)


# Assuming p.value is the column name for p-values in the sign section
###signif_asvs code NOT working -> is grouping non-signif ASVs... not sure why
View(inv_F_day1$sign)
inv_F_day1$sign$Feature_ID <- row.names(inv_F_day1$sign)
View(inv_F_day1$sign)

inv_F_spat$sign$Feature_ID <- row.names(inv_F_spat$sign)
View(inv_F_spat$sign)

# Extract the relevant columns from the 'sign' data frame
inv_F_sign_df_day1 <- data.frame(
  Feature_ID = inv_F_day1$sign$Feature_ID,
  p.value = inv_F_day1$sign$p.value,
  s.Control = inv_F_day1$sign$s.Control,
  s.High = inv_F_day1$sign$s.High,
  s.Low = inv_F_day1$sign$s.Low
)
inv_F_sign_df_spat <- data.frame(
  Feature_ID = inv_F_spat$sign$Feature_ID,
  p.value = inv_F_spat$sign$p.value,
  s.Control = inv_F_spat$sign$s.Control,
  s.High = inv_F_spat$sign$s.High,
  s.Low = inv_F_spat$sign$s.Low
)

# View the resulting data frame
View(inv_F_sign_df_day1)


#below create a "condition" column that will state which treatment is 
#enriched for an ASV for plots

# Initialize the "condition" column as an empty character vector
inv_F_sign_df_day1$condition <- ""
inv_F_sign_df_spat$condition <- ""

# Check each row for conditions and concatenate letters accordingly
inv_F_sign_df_day1$condition <- ifelse(inv_F_sign_df_day1$s.Control == 1, paste0(inv_F_sign_df_day1$condition, "C"), inv_F_sign_df_day1$condition)
inv_F_sign_df_day1$condition <- ifelse(inv_F_sign_df_day1$s.High == 1, paste0(inv_F_sign_df_day1$condition, "H"), inv_F_sign_df_day1$condition)
inv_F_sign_df_day1$condition <- ifelse(inv_F_sign_df_day1$s.Low == 1, paste0(inv_F_sign_df_day1$condition, "L"), inv_F_sign_df_day1$condition)

inv_F_sign_df_spat$condition <- ifelse(inv_F_sign_df_spat$s.Control == 1, paste0(inv_F_sign_df_spat$condition, "C"), inv_F_sign_df_spat$condition)
inv_F_sign_df_spat$condition <- ifelse(inv_F_sign_df_spat$s.High == 1, paste0(inv_F_sign_df_spat$condition, "H"), inv_F_sign_df_spat$condition)
inv_F_sign_df_spat$condition <- ifelse(inv_F_sign_df_spat$s.Low == 1, paste0(inv_F_sign_df_spat$condition, "L"), inv_F_sign_df_spat$condition)


#List signif ASVs ----

#Day 1 ASVs
# Filter rows where p.value is less than 0.05
significant_rows <- inv_F_sign_df_day1$sign[inv_F__df_day1$sign$p.value < 0.05, ]

# Extract the Feature_IDs from those rows
significant_ASVs_day1 <- significant_rows$Feature_ID

significant_ASVs_day1 <- significant_ASVs_day1[!is.na(significant_ASVs_day1)]

significant_ASVs_day1 <- as.list(significant_ASVs_day1)

#Spat ASVs
significant_rows <- inv_F_spat$sign[inv_F_spat$sign$p.value < 0.05, ]

significant_ASVs_spat <- significant_rows$Feature_ID

significant_ASVs_spat <- significant_ASVs_spat[!is.na(significant_ASVs_spat)]

significant_ASVs_spat <- as.list(significant_ASVs_spat)

# Find common ASVs between the time points
common_ASVs <- intersect(significant_ASVs_day1, significant_ASVs_spat)
#mb2021: ASV42, 188, 237, 240, 1164

#instead making manual list of signif (p<0.05) ASVs
#To extract as list of ASVs

#note: removed 319 (enriched in control and high sal)

# significant_asvs_indices_day1 <- c(119, 244, 192, 341, 510, 240, 428, 76, 450, 468,
#                                    484, 579, 385, 116, 199, 290, 54, 114, 60, 149, 339,
#                                    300, 200, 373, 326, 193, 659, 609, 318, 524, 102, 179, 
#                                    232, 38, 641, 150, 131, 92, 130, 257, 34, 499, 394)
# 
# 
# significant_asvs_indices_spat <- c(333, 507, 580, 494, 656, 371, 373, 623, 227, 283, 332, 571, 662, 385, 314, 380)

#sort ASVs from highest to lowest
# sorted_indices <- sort(significant_asvs_indices_day1, decreasing = FALSE)
# 
# sorted_indices <- sort(significant_asvs_indices_spat, decreasing = FALSE)
# 
# # Prepend "ASV" to each value in the list
# significant_asvs_names <- paste("ASV", sorted_indices, sep = "")
# print(significant_asvs_names)

#Get pseq data

pseq<- mb2021_filtered_NOT_rarefied
pseq <- subset_samples(pseq, !Age %in% c("3 dpf", "18 dpf", "Spat"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ps <- psmelt2(pseq) #long format
View(ps)


#Match condition column to ps ---- 

#make a new column in ps called condition and merge condition columns from 
#inv_F_sign_df to match up to FeatureID and populate condition column in ps

View(inv_F_sign_df)

match_result <- match(ps$FeatureID, inv_F_sign_df$Feature_ID)

#If have dataframe where only signif rows selected

match_result2 <- match(ps$FeatureID, inv_F_sign_df$selected_rows)

# Assign condition values from inv_F_sign_df to ps, handling NA values
ps$condition <- ifelse(is.na(match_result), NA, inv_F_sign_df$condition[match_result])

View(ps)

#formulate output dataframe - get signif asvs and log_abundance values
#arrange = to arrange in order highest to lowest ASV#

output <- ps %>%
  filter(FeatureID %in% significant_asvs_names) %>% 
  group_by(Treatment, FeatureID, Family.x, Order, Phylum, Class) %>%
  mutate(Log_Abundance = log(value))%>% 
  arrange(match(FeatureID, significant_asvs_names))
View(output)

str(output)

#Try plotting ASVs that are only signif more abundant in HS and LS to make plot cleaner

output <- output %>%
  filter(condition %in% c("H", "L"))

output <- output %>%
  filter(condition %in% c("C"))

#Colours ----
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

# Define your custom colors
my_colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", 
               "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
               "#ff33bb", "#ffcc00", "#00ffcc", "#0033ff", "#66ff33")

install.packages("randomcoloR")
library(randomcoloR)
n <- 50
palette <- distinctColorPalette(n)


custom_shapes <- c(16, 15, 18, 17, 19, 20)  # Choose from a list of available shapes (0-25) in ggplot2

#order x-axis from highest to lowest ASVs
output$FeatureID <- factor(output$FeatureID, levels = unique(output$FeatureID))


#scatterplot signif ASVs ----
#eventually, want to add phylogenetic tree on x axis

ggplot(output, aes(x = FeatureID, y = Log_Abundance, color = Family.x, shape = condition)) + 
  geom_point(size = 6) +
  theme_classic() +
  facet_grid(~Treatment) +
  scale_color_manual(values=palette) +
  scale_shape_manual(values = c("C" = 17, "H" = 17, "L" = 15, "CL" = 6, "HL" = 1)) +  # Custom shapes for each condition
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.background = element_rect(color = "black", fill = NA),  # Add border outline
    legend.box.background = element_rect(color = NA),  # Remove legend border
    legend.key = element_rect(color = NA),  # Remove border around legend values
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 14)  # Adjust facet label font size
  ) + labs(y = "Log Relative Abundance")


View(output)

ggplot(output, aes(x = FeatureID, y = Log_Abundance, color = Order, shape = Treatment)) + 
  geom_point(size = 3) +
  theme_bw() +
  #facet_grid(~Treatment) +
  scale_color_manual(values=mycolors) +
  #scale_shape_manual(values = c("C" = 6, "H" = 2, "L" = 6, "CH" = 4, "CL" = 6, "HL" = 1)) +  # Custom shapes for each condition
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size=7, angle=45, hjust=1)
  )


ggplot(output, aes(x = FeatureID, y = Log_Abundance, color = Order, shape = Phylum)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=col_vector) +
  facet_grid(~Treatment) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_shape_manual(values = custom_shapes)

ggplot(output, aes(x = FeatureID, y = Log_Abundance, color = Order, shape = Phylum)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values=col_vector) +
  facet_grid(~Treatment) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_shape_manual(values = custom_shapes)




#To make ASV significance as a matrix for heatmap 
significant_asvs_matrix <- inv_F$str[significant_asvs_indices, ]
View(significant_asvs_matrix)
str(significant_asvs_indices)

pheatmap(significant_asvs_matrix, angle_col="45")

#exclude certain columns
selected_columns <- c("Day 01", "Day 03", "Day 06", "Day 15", "Spat")  # Replace with your actual column names
selected_columns <- c("Day 01", "Day 03", "Day 06", "Day 08", "Day 15")

selected_columns <- c("Control", "High salinity", "Low salinity")

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


