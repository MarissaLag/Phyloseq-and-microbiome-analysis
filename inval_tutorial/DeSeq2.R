#DESeq2 test for differential species abundance test
#Code from Elliot Scanes


#Should I be normalizing the data beforehand???

library(phyloseq)
library(microbiome)

#load phyloseq object
pseq <- Marissa_mb2021_filtered_20240203

pseq <- mb2021_filtered_NOT_rarefied

pseq <- MU42022_filtered_NOT_rarefied

pseq <- MU42022_filtered_NOT_rarefied_moreTAXA 

pseq <- Marissa_MU42022_rare_nochloro

  
#filter samples

#MU42022
pseq <- subset_samples(pseq, !Genetics %in% "4")
  
pseq <- subset_samples(pseq, Age %in% "1 dpf")
pseq <- subset_samples(pseq, !Library_Name %in% c("T9r1", "T9r3"))

pseq <- subset_samples(pseq, Age %in% "Day 15")
pseq <- subset_samples(pseq, !Family %in% "9")
pseq <- subset_samples(pseq, !Organism %in% "Algae")
View(pseq@sam_data)


# Agglomerating family names (mb2021)
pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)]  <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4
pseq@sam_data$Family <- as.factor(pseq@sam_data$Family)


#Core

pseq <- core(pseq, detection = .1/100, prevalence = 90/100)

View(pseq@otu_table)

 #Install deseq
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager") #download deseq
# 
# BiocManager::install("DESeq2") #Install deseq
# 
# BiocManager::install("DESeq2", force = TRUE)
# 
# BiocManager::install("GenomeInfoDb")
# BiocManager::install("DESeq2")

#filter out probiotic (ASV7 and 18) for MU42022
taxa_to_remove <- c("ASV7", "ASV18")

# Create a logical vector indicating which taxa to keep
taxa_to_keep <- !(taxa_names(pseq) %in% taxa_to_remove)

# Prune the taxa from the phyloseq object
pseq <- prune_taxa(taxa_to_keep, pseq)

library(DESeq2) #load deseq

#Should relative abundance be used? Get error when I try

#pseq_rel <- microbiome::transform(pseq, "compositional")


DeSeq <- phyloseq_to_deseq2(pseq, ~ Treatment) #convert phyloseq to deseq object

DeSeq2 <- DESeq(DeSeq) #run deseq analysis

#Notes about the DESeq fxn:
#This function performs a default analysis through the steps:
#estimation of size factors: estimateSizeFactors #if data is rarefied this will be null (all same size)
#estimation of dispersion: estimateDispersions
#Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest #abundance data normalization


#format for code below must be: contrast = c('factorName', 'numeratorLevel', 'denominatorLevel')

res_high_vs_control <- results(DeSeq2, contrast = c("Treatment", "High salinity", "Control"))
res_low_vs_control <- results(DeSeq2, contrast = c("Treatment", "Low salinity", "Control"))
#res_high_vs_low <- results(DeSeq2, contrast = c("Treatment", "High salinity", "Low salinity"))

#Family effects
# res_F1_vs_F2 <- results(DeSeq2, contrast = c("Family", "1", "2"))
# res_F2_vs_F3 <- results(DeSeq2, contrast = c("Family", "2", "3"))
# res_F3_vs_F4 <- results(DeSeq2, contrast = c("Family", "3", "4"))
# res_F1_vs_F3 <- results(DeSeq2, contrast = c("Family", "1", "3"))
# res_F2_vs_F4 <- results(DeSeq2, contrast = c("Family", "2", "4"))
# res_F1_vs_F4 <- results(DeSeq2, contrast = c("Family", "1", "4"))



res_PB_vs_control <- results(DeSeq2, contrast = c("Treatment", "Probiotics", "Control"))
res_PBH_vs_control <- results(DeSeq2, contrast = c("Treatment", "Probiotics + HT", "Control"))
res_HT_vs_control <- results(DeSeq2, contrast = c("Treatment", "High temperature", "Control"))

** from this point it’s optional but this code helps filter the results **

res_dat_low <- cbind(as(res_low_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_low_vs_control), ], "matrix")) #make the results a data frame
res_dat_high <- cbind(as(res_high_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_high_vs_control), ], "matrix")) #make the results a data frame
#res_dat_high_low <- cbind(as(res_high_vs_low, "data.frame"), as(tax_table(pseq)[rownames(res_high_vs_low), ], "matrix")) #make the results a data frame


res_dat_PB <- cbind(as(res_PB_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_PB_vs_control), ], "matrix")) #make the results a data frame
res_dat_PBH <- cbind(as(res_PBH_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_PBH_vs_control), ], "matrix")) #make the results a data frame
res_dat_HT <- cbind(as(res_HT_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_HT_vs_control), ], "matrix")) #make the results a data frame


alpha = 0.05 #Set alpha

sigtab_high = res_dat_high[which(res_dat_high$padj < alpha), ] #filter out significant results
sigtab_low = res_dat_low[which(res_dat_low$padj < alpha), ]
sigtab_high_low = res_dat_high_low[which(res_dat_high_low$padj < alpha), ]

sigtab_PB = res_dat_PB[which(res_dat_PB$padj < alpha), ] #filter out significant results
sigtab_PBH = res_dat_PBH[which(res_dat_PBH$padj < alpha), ]
sigtab_HT = res_dat_HT[which(res_dat_HT$padj < alpha), ]

#print only the significant results                              
sigtab_high
sigtab_low
#sigtab_high_low


#plotting
#source: https://joey711.github.io/phyloseq-extensions/DESeq2.html

library("ggplot2")
library(viridis)
library(hrbrthemes)
theme_set(theme_bw())


# Phylum order
x = tapply(sigtab_high$log2FoldChange, sigtab_high$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab_high$Order = factor(as.character(sigtab_high$Class), levels=names(x))
# Genus order
x = tapply(sigtab_high$log2FoldChange, sigtab_high$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_high$Genus = factor(as.character(sigtab_high$Genus), levels=names(x))

mycolors <- c("#E69F00", "#CC79A7", "#009E73", "#56B4E9", "#F0E442", "#999999")


ggplot(sigtab_PB, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=7) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold")) +
  labs(title = "High salinity vs. Control", x = "", y = "Log2-Fold-Change") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = mycolors)

  scale_color_manual(values = c("Alteromonadaceae" = "pink", "Halieaceae" = "#377EB8", "Rhodobacteraceae" = "#4DAF4A", "Saprospiraceae" = "yellow"))

ggplot(sigtab_high, aes(x=Family, y=log2FoldChange, color=Family)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) +
  labs(title = "High salinity vs. Control", x = "", y = "Log2-Fold-Change") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid = element_blank()) +
  scale_color_ipsum()


##Horizontal var graph -----


#Signif ASVs between High sal and control treatment
View(sigtab_high)

#Signif ASVs between Low sal and control treatment
View(sigtab_low)

#Need to convert row names (ASV#s) into a column)
library(tibble)
library(RColorBrewer)

# Convert row names to a column named "ASV"
sigtab_low <- rownames_to_column(sigtab_low, var = "ASV")
sigtab_high <- rownames_to_column(sigtab_high, var = "ASV")

sigtab_PB <- rownames_to_column(sigtab_PB, var = "ASV")
sigtab_PBH <- rownames_to_column(sigtab_PBH, var = "ASV")
sigtab_HT <- rownames_to_column(sigtab_HT, var = "ASV")

#Want family identity info on bar graph too (but colouring is too much) so making a new column 
sigtab_high$Combined_Info <- paste(sigtab_high$Genus, sigtab_high$ASV, sep = ";")
sigtab_low$Combined_Info <- paste(sigtab_low$Genus, sigtab_low$ASV, sep = ";")

sigtab_PB$Combined_Info <- paste(sigtab_PB$Genus, sigtab_PB$ASV, sep = ";")
sigtab_PBH$Combined_Info <- paste(sigtab_PBH$Genus, sigtab_PBH$ASV, sep = ";")
sigtab_HT$Combined_Info <- paste(sigtab_HT$Genus, sigtab_HT$ASV, sep = ";")

custom_palette <- brewer.pal(12, "Set3")

custom_palette_low <- c(  "#FFFFB3", "#6a3d9a", "#00ffcc","green", "#BEBADA", "brown", "#80B1D3","#FDB462","#CCEBC5","red","#FFED6F", "#B3DE69","#FB8072", "#80B1D3",   
                    "yellow", "#8DD3C7", "#FCCDE5", 
                    "red", "#BC80BD", "#FFED6F", "#ff7f00", "#a6cee3","#ffcc00" )

custom_palette <- c(  "#FFFFB3", "#00ffcc","#BEBADA", "brown", "#B3DE69","#FB8072", "#80B1D3","#FDB462","#CCEBC5",   
                      "yellow", "#8DD3C7", "#FCCDE5", 
                      "red", "#BC80BD", "#FFED6F", "#ff7f00", "#a6cee3","#ffcc00" )



custom_palette <- c("#a6cee3", "#b2df8a","#ff7f00","#fb9a99",  "#fdbf6f","#e31a1c" ,"lightyellow","#6a3d9a", 
              "#fdbf6f", "#cab2d6","brown", "#ffcc00", "#00ffcc", "grey", "skyblue3", "#66ff33" )


custom_palette <- c("#33a02c", "#a6cee3","#ffcc00", "brown", "#cab2d6", "#ff7f00", "#b2df8a", "#e31a1c","#6a3d9a", 
                     "#fb9a99","#ff33bb", "#00ffcc", "#66ff33", "#a6cee3", "#b2df8a","#66ff33","#ff7f00","#fb9a99",  "#fdbf6f", "#ff7f00", "#e31a1c","#6a3d9a", 
                    "#fdbf6f", "#cab2d6","#ff33bb", "#ffcc00", "#00ffcc")



ggplot(sigtab_low, aes(x = reorder(Combined_Info, log2FoldChange), y = log2FoldChange, fill = Family)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  labs(title = "", 
       x = "Genus; ASV", 
       y = "Log2 Fold Change") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),   # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = custom_palette_low)

#Filter out ASVs shared between HS and LS

# Identify the shared ASVs
shared_asvs <- intersect(sigtab_high$ASV, sigtab_low$ASV)

# Filter out the shared ASVs from both datasets
sigtab_high_filtered <- sigtab_high %>%
  filter(!ASV %in% shared_asvs)

sigtab_low_filtered <- sigtab_low %>%
  filter(!ASV %in% shared_asvs)

head(sigtab_high_filtered)
head(sigtab_low_filtered)

ggplot(sigtab_low_filtered, aes(x = reorder(Combined_Info, log2FoldChange), y = log2FoldChange, fill = Class)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  labs(title = "1 dpf - Significant ASVs - Low salinity vs Control", 
       x = "Family; ASV", 
       y = "Log2 Fold Change") +
  theme(panel.grid.major = element_blank(),   # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = custom_palette)


#Venn diagram of shared ASVs ----

# Extract ASV names from each set
asv_PB <- sigtab_PB$ASV
asv_PBH <- sigtab_PBH$ASV
asv_HT <- rownames(sigtab_HT)

# Create a list of ASV sets
asv_list <- list(
  PB = asv_PB,
  PBH = asv_PBH,
  HT = asv_HT
)

mycols <- c(PB = "lightgreen", PBH = "lightblue", HT = "#F8766D")

# Create Venn diagram
venn_result <- venn.diagram(
  x = list(
    "Probiotics" = asv_list[["PB"]],
    "Probiotics + HT" = asv_list[["PBH"]],
    "High temperature (HT)" = asv_list[["HT"]]
  ),
  category.names = c("Probiotics", "Probiotics + HT", "High temperature (HT)"),
  fill = mycols,
  filename = NULL
)

# Plot the Venn diagram
grid.draw(venn_result)








#Graph certain ASVs ----

library(ggplot2)
library(dplyr)

psmelt <- psmelt(pseq)

# Assuming psmelt is already loaded in your environment
# Filter for ASV40
Alii_data <- psmelt %>%
  filter(Genus == "Aliikangiella")

# Create the box plot
ggplot(Alii_data, aes(x = Treatment, y = Abundance)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Average Abundance of ASV40 by Treatment",
       x = "Treatment",
       y = "Abundance") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  facet_wrap(~OTU)

#Normalize data for figures ---- 
#Data normalisation using DESeq2:
#counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") #download deseq

BiocManager::install("DESeq2") #Install deseq

BiocManager::install("DESeq2", force = TRUE)

BiocManager::install("GenomeInfoDb")
BiocManager::install("DESeq2")

library(DESeq2)
library(phyloseq)


pseq <- MU42022_filtered_NOT_rarefied
pseq <- mb2021_filtered_NOT_rarefied
#how to pick which factors to use?? Checked, does not matter which factors you use, will normalize the same. 
#need to remove NAs (algae)
DeSeq <- phyloseq_to_deseq2(pseq, ~ Age) #convert phyloseq to deseq object

DeSeq2 <- DESeq(DeSeq)

fact1 = sample_data(pseq)
fact = as.matrix.data.frame(fact1)
fact = as.data.frame(fact)

count1 = otu_table(pseq)
count = as.matrix.data.frame(count1)
count = as.data.frame(count)

# # two-way analysis with 2 levels each (treatment + time)
# # countData = ASV counts with samples in columns  
# # colData = Sample meta data (factors, etc)  
# # design = the right hand side of a GLM formula

ASV_deseq <- DESeqDataSetFromMatrix(countData = t(count),
                                    colData = fact,
                                    design = ~ Age)

# positive counts normalisaton (accounts for zero-inflation)
ASV_deseq <- estimateSizeFactors(DeSeq2, type = "poscounts")
sizeFactors(ASV_deseq)

# size-factor corrected data are calculated by dividing the raw counts by the sample size factor and adding 0.5 to correct for the zeros
ASV_deseq_norm <- sapply(row.names(ASV_deseq), function(x){
  plotCounts(ASV_deseq, x, "Age", returnData = TRUE, normalized = TRUE)$count
  
})
rownames(ASV_deseq_norm) <- colnames(ASV_deseq)

# remove the addition of 0.5 to all entries
ASV_deseq_norm <- ASV_deseq_norm - 0.5
ASV_deseq_norm[1:10, 1:10]

# make dataset with integers
ASV_deseq_norm <- round(ASV_deseq_norm)
ASV_deseq_norm[1:10, 1:10]

View(ASV_deseq_norm)

# check for zeros across all samples (produced by normalisation)
dim(ASV_deseq_norm[, colSums(ASV_deseq_norm) == 0]) # 0

# check singletons
dim(ASV_deseq_norm[, colSums(ASV_deseq_norm) == 1])

#save normalised table
write.csv(ASV_deseq_norm, "~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files//normalised_ASV_table_mb2021_Age.csv")

#Make new pseq object with normalized ASV table
mb2021_filtered_NOT_rarefied <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/mb2021_filtered_NOT_rarefied.rds")

#Check format
View(pseq@otu_table)
View(normalised_ASV_table_mb2021_Age)

#correct to make forst column as row names

normalised_ASV_table_mb2021_Age <- as.data.frame(normalised_ASV_table_mb2021_Age)

rownames(normalised_ASV_table_mb2021_Age) <- normalised_ASV_table_mb2021_Age[, 1]

# Remove the first column from the data frame
normalised_ASV_table_MU42022 <- normalised_ASV_table_MU42022[, -1]

# Create a new otu_table object
new_otu_table <- otu_table(normalised_ASV_table_MU42022, taxa_are_rows = FALSE)
View(new_otu_table)
str(new_otu_table)
# Replace the old OTU table with the new one
#otu_table(mb2021_filtered_NOT_rarefied) <- new_otu_table
otu_table(MU42022_filtered_NOT_rarefied) <- new_otu_table

#If giving error that sample names do not match
# Extract the sample names from new_otu_table
new_sample_names <- sample_names(new_otu_table)

# Verify the replacement
View(mb2021_filtered_NOT_rarefied@otu_table)
#Save as different RDS file
saveRDS(mb2021_filtered_NOT_rarefied, file = "mb2021_filtered_NOT_rarefied_normalized.rds")



