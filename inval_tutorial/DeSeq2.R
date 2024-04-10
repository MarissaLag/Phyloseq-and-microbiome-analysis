#DESeq2 test for differential species abundance test
#Code from Elliot Scanes


#Should I be normalizing the data beforehand???


#load phyloseq object
pseq <- Marissa_mb2021_filtered_20240203

#filter samples
pseq <- subset_samples(pseq, !Age %in% "3 dpf")
pseq <- subset_samples(pseq, !Library_Name %in% c("T9r1", "T9r3"))

#Select time-point
pseq <- subset_samples(pseq, Age %in% "Spat")
View(pseq@sam_data)

#Can I look at core member? If so,

pseq<- core(pseq, detection = .1/100, prevalence = 90/100)

View(pseq@otu_table)

 #Install deseq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") #download deseq

BiocManager::install("DESeq2") #Install deseq

BiocManager::install("DESeq2", force = TRUE)

BiocManager::install("GenomeInfoDb")
BiocManager::install("DESeq2")


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
res_high_vs_low <- results(DeSeq2, contrast = c("Treatment", "High salinity", "Low salinity"))


** from this point it’s optional but this code helps filter the results **

res_dat_low <- cbind(as(res_low_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_low_vs_control), ], "matrix")) #make the results a data frame
res_dat_high <- cbind(as(res_high_vs_control, "data.frame"), as(tax_table(pseq)[rownames(res_high_vs_control), ], "matrix")) #make the results a data frame
res_dat_high_low <- cbind(as(res_high_vs_low, "data.frame"), as(tax_table(pseq)[rownames(res_high_vs_low), ], "matrix")) #make the results a data frame


alpha = 0.05 #Set alpha

sigtab_high = res_dat_high[which(res_dat_high$padj < alpha), ] #filter out significant results
sigtab_low = res_dat_low[which(res_dat_low$padj < alpha), ]
sigtab_high_low = res_dat_high_low[which(res_dat_high_low$padj < alpha), ]

#code below not working but taxa already included?
#add the taxonomy back in
sigtab_high = cbind(as(sigtab_high, "data.frame"), as(tax_table(pseq[rownames(sigtab_high), ], "matrix")) 
sigtab_low = cbind(as(sigtab_low, "data.frame"), as(tax_table(pseq[rownames(sigtab_low), ], "matrix")) #add the taxonomy back in

#print only the significant results                              
sigtab_high
sigtab_low
sigtab_high_low


#plotting
#source: https://joey711.github.io/phyloseq-extensions/DESeq2.html

library("ggplot2")
library(viridis)
library(hrbrthemes)
theme_set(theme_bw())


# Phylum order
x = tapply(sigtab_high$log2FoldChange, sigtab_high$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_high$log2FoldChange, sigtab_high$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_high$Genus = factor(as.character(sigtab_high$Genus), levels=names(x))

mycolors <- c("#E69F00", "#CC79A7", "#009E73", "#56B4E9", "#F0E442", "#999999")
mycolors <- c("#F0E442")


ggplot(sigtab_high_low, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=7) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold")) +
  labs(title = "High salinity vs. Low salinity", x = "", y = "Log2-Fold-Change") +
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
