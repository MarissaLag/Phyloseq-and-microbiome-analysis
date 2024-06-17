#DESeq2 test for differential species abundance test
#Code from Elliot Scanes


#Should I be normalizing the data beforehand???


#load phyloseq object
pseq <- Marissa_mb2021_filtered_20240203

pseq <- mb2021_filtered_NOT_rarefied

pseq <- MU42022_filtered_NOT_rarefied

pseq <- MU42022_filtered_NOT_rarefied_moreTAXA 

pseq <- Marissa_MU42022_rare

pseq <- mb2021_filtered_NOT_rarefied_moreTAXA_5percent_removal
  
#filter samples

#MU42022
pseq <- subset_samples(pseq, !Genetics %in% "4")
  
pseq <- subset_samples(pseq, !Age %in% "3 dpf")
pseq <- subset_samples(pseq, !Library_Name %in% c("T9r1", "T9r3"))

#Select time-point
pseq <- subset_samples(pseq, Age %in% "1 dpf")
pseq <- subset_samples(pseq, !Family %in% "9")
pseq <- subset_samples(pseq, !Organism %in% "Algae")
View(pseq@sam_data)


# Agglomerating family names (mb2021)
pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)]  <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4
pseq@sam_data$Family <- as.factor(pseq@sam_data$Family)


#Can I look at core member? If so,

pseq<- core(pseq, detection = .1/100, prevalence = 90/100)

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
#res_HT_vs_control <- results(DeSeq2, contrast = c("Treatment", "High temperature", "Control"))

** from this point it’s optional but this code helps filter the results **

res_dat_low <- cbind(as(res_low_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_low_vs_control), ], "matrix")) #make the results a data frame
res_dat_high <- cbind(as(res_high_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_high_vs_control), ], "matrix")) #make the results a data frame
#res_dat_high_low <- cbind(as(res_high_vs_low, "data.frame"), as(tax_table(pseq)[rownames(res_high_vs_low), ], "matrix")) #make the results a data frame


res_dat_PB <- cbind(as(res_PB_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_PB_vs_control), ], "matrix")) #make the results a data frame
res_dat_PBH <- cbind(as(res_PBH_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_PBH_vs_control), ], "matrix")) #make the results a data frame
#res_dat_HT <- cbind(as(res_HT_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_HT_vs_control), ], "matrix")) #make the results a data frame


alpha = 0.05 #Set alpha

sigtab_high = res_dat_high[which(res_dat_high$padj < alpha), ] #filter out significant results
sigtab_low = res_dat_low[which(res_dat_low$padj < alpha), ]
#sigtab_high_low = res_dat_high_low[which(res_dat_high_low$padj < alpha), ]

sigtab_PB = res_dat_PB[which(res_dat_PB$padj < alpha), ] #filter out significant results
sigtab_PBH = res_dat_PBH[which(res_dat_PBH$padj < alpha), ]
#sigtab_HT = res_dat_HT[which(res_dat_HT$padj < alpha), ]

#code below not working but taxa already included?
#add the taxonomy back in
sigtab_high = cbind(as(sigtab_high, "data.frame"), as(tax_table(pseq[rownames(sigtab_high), ], "matrix")) 
sigtab_low = cbind(as(sigtab_low, "data.frame"), as(tax_table(pseq[rownames(sigtab_low), ], "matrix")) 

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


ggplot(sigtab_high, aes(x=Genus, y=log2FoldChange, color=Class)) + geom_point(size=7) + 
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

#Want family identity info on bar graph too (but colouring is too much) so making a new column 
sigtab_high$Combined_Info <- paste(sigtab_high$Family, sigtab_high$ASV, sep = ";")

sigtab_low$Combined_Info <- paste(sigtab_low$Family, sigtab_low$ASV, sep = ";")

custom_palette <- brewer.pal(12, "Set3")

custom_palette <- c( "#FFFFB3", "#8DD3C7", "#FB8072", "#FDB462",
                    "#80B1D3", "#B3DE69", "#FCCDE5", 
                    "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#BEBADA")


custom_palette <- c(  "#a6cee3", "#b2df8a", "#fdbf6f", "#ff7f00", "#e31a1c","#6a3d9a", 
             "#fb9a99", "#fdbf6f", "#cab2d6","#ff33bb", "#ffcc00", "#00ffcc", "#66ff33")


custom_palette <- c("#33a02c", "#a6cee3","#ffcc00", "brown", "#cab2d6", "#ff7f00", "#b2df8a", "#e31a1c","#6a3d9a", 
                     "#fb9a99","#ff33bb", "#00ffcc", "#66ff33")



ggplot(sigtab_low, aes(x = reorder(Combined_Info, log2FoldChange), y = log2FoldChange, fill = Class)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  labs(title = "Spat - Significant ASVs - High salinity vs Control", 
       x = "Genus; ASV", 
       y = "Log2 Fold Change") +
  theme(panel.grid.major = element_blank(),   # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = custom_palette)

