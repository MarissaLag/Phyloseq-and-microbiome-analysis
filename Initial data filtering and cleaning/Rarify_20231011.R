##testing MDJ.rds data


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

library("devtools")
library(phyloseq)
library(microbiome)

pseq<- readRDS("Marissa_.rds")

#pseq <- Marissa_Oyster

pseq <- Marissa_MU42022_rare_nochloro

pseq

##Want to seperate Marissa's and James' samples

excluded_samples <-

Rare_filtered <- Rare[, !colnames(Rare) %in% excluded_samples]



#create objects

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data

# check if any OTUs are not present in any samples
any(taxa_sums(pseq) == 0)

#source for removing chloro/mito/archaea; https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html#making-a-phyloseq-object

#remove chloroplast, mito, chimera

pseq1 <- subset_taxa(pseq,Class!="c__Chloroplast")

pseq2 <- subset_taxa(pseq1,Order!="o__Mitochondria")

ps1 <- subset_taxa(pseq2,Kingdom!="k__Archaea")

#chloroplast still in data under order level - remove?

pseq1 <- subset_taxa(pseq,Order!="Chloroplast")

ps1 <- pseq1

rank_names(ps1)



###remove low prev

plot(sort(taxa_sums(x1), TRUE), type="h", ylim=c(0, 10000))

x1 = prune_taxa(taxa_sums(ps1) > 200, ps1) 
x2 = prune_taxa(taxa_sums(ps1) > 500, ps1) 
x3 = prune_taxa(taxa_sums(ps1) > 1000, ps1)

##using x2

library(microbiome)
summarize_phyloseq(ps1)
summarize_phyloseq(x1)


*************************************************
#rarify - removing abundance code with below code does not seem to work...not sure why. Ignore for now.

sequence_abundance <- taxa_sums(ps1)

threshold <- 0.005

# Get the names of sequences that meet the threshold
abundant_sequences <- names(sequence_abundance[sequence_abundance > threshold])

# Filter out sequences not meeting the threshold and create a new phyloseq object
physeq_rarified


keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps1)

*************************************************

##rarify data to make sequences an even read depth - selecting read depth of 10,000 = any samples with fewer than 10,000 total reads will be removed, all samples will be equalized to 5000 reads for Denman's samples, 10,000 for Marissa/James'

Rare <-rarefy_even_depth(x2, sample.size= 5000)

#Marissa/James = 17 samples removed because they contained fewer reads than `sample.size' - first 5 reads are T13-1,T13-2-2,T14-2-2,T15-S-1,T16-10-r3

##check if seq depth is 10,000 or 5,000

sample_depths <- sample_sums(pseq)

print(sample_depths)

##yes, all samples have 10,000 or 5,000 reads now.

##rename Rare to ps1

ps1 <- Rare

###calculate relative abundance at family, phylum, and genus level.

pseq_fam <- microbiome::aggregate_rare(ps1, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(ps1, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(ps1, level = "Genus", detection = 50/100, prevalence = 70/100)

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")



##Get a quick look at the data (plot ordination)

set.seed(4235421)

ord <- ordinate(ps1, "MDS", "bray")

#plot MDS/PcoA - can set "colour" and "shape" for any of your variableas
#geompoint controls data point size on plot

plot_ordination(ps1, ord, color = "Family.1", shape = "Percent.Fouled") + geom_point(size = 4) + scale_shape_binned()

plot_ordination(ps1, ord, color = "Factor") + geom_point(size = 4)

plot_ordination(ps1, ord, color = "Family.1") + geom_point(size = 4)

plot_ordination(ps1, ord, color = "Treatment", shape = "Age") + geom_point(size = 4)


p <- plot_ordination(ps1, ord, color = "Factor") + geom_point(size = 4)

