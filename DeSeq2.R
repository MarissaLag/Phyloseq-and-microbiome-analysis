#deSeq2 differential abundance analysis
source: https://sw1.github.io/teaching/phyloseq.html
source: https://joey711.github.io/phyloseq-extensions/DESeq2.html


#Many issues with this code that needs to be fixed
#attempting to run Deseq2 differential abundance analysis and/or ANCOM


Marissa_mb2021_filtered_20240203 <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Marissa_mb2021_filtered_20240203.rds")

PS <- Marissa_mb2021_filtered_20240203

PS <- subset_samples(PS, !Age %in% c("3 dpf", "Spat", "18 dpf"))


install.packages("DESeq2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(phyloseq)
library(tidyverse)
library(ape)
library(DESeq2)

PS3 <- prune_samples(!is.na(sample_data(PS)$Treatment),PS)
PS3 <- filter_taxa(PS3,function(x) sum(x) > 0,prune = TRUE)

diagdds <- phyloseq_to_deseq2(PS3, ~ Treatment)


## converting counts to integer mode
diagdds <- DESeq(diagdds, test='Wald', fitType='parametric')

res <- results(diagdds, cooksCutoff = FALSE)
res

#only compared Age Spat vs 1.dpf I think...


res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PS3)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)


#visualization

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


#running another tutorial
#source: https://www.yanh.org/2021/01/01/microbiome-r/#differential-abundance-analysis
#source: https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html#deseq2


ps <- Marissa_mb2021_filtered_20240203 
 

sample_data(ps)$Treatment <- as.factor(sample_data(ps)$Treatment) # factorize for DESeq2

ps.taxa <- tax_glom(ps, taxrank = 'Genus', NArm = FALSE)

# pairwise comparison between Treatments
ps.taxa.sub <- subset_samples(ps.taxa, Treatment %in% c("Control", "High salinity", "Low salinity"))

# filter sparse features, with > 90% zeros ##NOT DONE: code below not working

ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

ps_ds = phyloseq_to_deseq2(ps.taxa.pse.sub, ~ body.site)




# use alternative estimator on a condition of "every gene contains a sample with a zero"
ds <- estimateSizeFactors(ps_ds, type="poscounts")
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05 
res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
taxa_sig = rownames(res[1:20, ]) # select bottom 20 with lowest p.adj values
ps.taxa.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
# Only keep gut and tongue samples
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa.pse.sub)), ps.taxa.rel.sig)



# Agglomerates data to Genus level
load.packages(mia)


tse_genus <- agglomerateByRank(ps, rank = "Genus")
tse_genus <- tax_glom(ps, taxrank = 'Genus', NArm = FALSE)

# Perform clr transformation. A Pseudocount of 1 needs to be added, 
# because the data contains zeros and the clr transformation includes a 
# log transformation.
tse_genus <- transformAssay(tse_genus, method = "clr", pseudocount = 1)

ds2 <- DESeqDataSet(tse_genus, ~Treatment)

