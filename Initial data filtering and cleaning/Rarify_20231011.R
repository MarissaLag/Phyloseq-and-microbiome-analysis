##testing MDJ.rds data

#Load packages ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

library("devtools")
library(phyloseq)
library(microbiome)

#Load un-rarefied data ----

pseq<- readRDS("Marissa_.rds")

pseq <- Marissa_mb2021_filtered_20240203

pseq <- Marissa_mb2021_unfiltered

pseq

#removing day 3 right away

pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))

#create objects

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data

#Sanity check
#check if any OTUs are not present in any samples (want false)
any(taxa_sums(pseq) == 0)

#if true

pseq_filtered <- prune_taxa(taxa_sums(pseq) > 0, pseq)
any(taxa_sums(pseq_filtered) == 0)

pseq <- pseq_filtered

#source for removing chloro/mito/archaea; https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html#making-a-phyloseq-object

#Remove chloroplast, mito, archaea, chimera ----

pseq1 <- subset_taxa(pseq,Class!="c__Chloroplast")

pseq2 <- subset_taxa(pseq1,Order!="o__Mitochondria")

ps1 <- subset_taxa(pseq2,Kingdom!="k__Archaea")

#chloroplast still in data under order level - remove?

pseq1 <- subset_taxa(pseq,Order!="Chloroplast")

ps1 <- pseq1

rank_names(ps1)

#Remove low prev ----

plot(sort(taxa_sums(x2), TRUE), type="h", ylim=c(0, 10000))

x1 = prune_taxa(taxa_sums(ps1) > 200, ps1) 
x2 = prune_taxa(taxa_sums(ps1) > 500, ps1) 
x3 = prune_taxa(taxa_sums(ps1) > 1000, ps1)

##using x2

library(microbiome)
summarize_phyloseq(ps1)
summarize_phyloseq(x2)


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

#Depending on data, may not always have to rarefy read counts
#if rarefication curves look good, don't rarefy. If there are one or two outliers, remove them and don't rarefy

#Also Jadi normalizes her data as well (even if rarefied or not) with DESEq2 negative binomial dist'd model

Marissa_MU42022_unfiltered <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/Marissa_MU42022_unfiltered.rds")

Marissa_mb2021_unfiltered <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/Marissa_mb2021_unfiltered.rds")

pseq <- Marissa_mb2021_unfiltered

pdf("Phyloseq and microbiome analysis/Old RDA files/rarefaction_curve.pdf")
raref.curve <- rarecurve(ASV_table_data, ylab = "ASV Count") # Long step
dev.off()

#Check if data needs rarefying
#Rarefication curves ----

pseq <- x2

#check reads

library(ggplot2) 
library(data.table)

# Check read count
readcount = data.table(as(sample_data(pseq), "data.frame"),
                         TotalReads = sample_sums(pseq), 
                         keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID")

ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")

head(readcount[order(readcount$TotalReads), c("SampleID", "TotalReads")])

#mb2021: Remove samples F4L18, T10r3, T9r2 (did not work at all)

pseq <- subset_samples(pseq, !Sample.ID %in% c("F4L18", "T10"))

otu.rare = otu_table(pseq)
otu.rare = as.data.frame(otu.rare)
sample_names = rownames(otu.rare)

#Generate rarefaction curve, rarefaction curve could be used to determined whether the sequencing depth cover microbial diversity of the sample.
# we will use vegan rarecurve 
library(vegan)

otu.rarecurve <- rarecurve(otu.rare, step = 10000, label = FALSE, xlim = c(0, 100000))

raref.curve <- rarecurve(otu.rare, label = FALSE, ylab = "ASV Count")

##rarify data to make sequences an even read depth - selecting read depth of 10,000 = any samples with fewer than 10,000 total reads will be removed, all samples will be equalized to 5000 reads for Denman's samples, 10,000 for Marissa/James'

Rare <-rarefy_even_depth(x2, sample.size= 3000)

#Marissa/James = 17 samples removed because they contained fewer reads than `sample.size' - first 5 reads are T13-1,T13-2-2,T14-2-2,T15-S-1,T16-10-r3

##check if seq depth is 10,000 or 5,000

#mb2021 rarefied to 3000

sample_depths <- sample_sums(pseq)

print(sample_depths)

##yes, all samples have 10,000 or 5,000 reads now.

##rename Rare to ps1

ps1 <- Rare

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq, "compositional")

pseq.phy.rel <- microbiome::transform(pseq, "compositional")

pseq.gen.rel <- microbiome::transform(pseq, "compositional")

pseq.core <- microbiome::transform(pseq.core, "compositional")



pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq_fam <- microbiome::aggregate_rare(ps1, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(ps1, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(ps1, level = "Genus", detection = 50/100, prevalence = 70/100)

##Get a quick look at the data (plot ordination)

set.seed(4235421)

ord <- ordinate(pseq, "MDS", "bray")

#plot MDS/PcoA - can set "colour" and "shape" for any of your variableas
#geompoint controls data point size on plot

plot_ordination(ps1, ord, color = "Family.1", shape = "Percent.Fouled") + geom_point(size = 4) + scale_shape_binned()

plot_ordination(ps1, ord, color = "Factor") + geom_point(size = 4)

plot_ordination(ps1, ord, color = "Family.1") + geom_point(size = 4)

plot_ordination(pseq, ord, color = "Treatment", shape = "Age") + geom_point(size = 4)


p <- plot_ordination(pseq, ord, color = "Factor") + geom_point(size = 4)

