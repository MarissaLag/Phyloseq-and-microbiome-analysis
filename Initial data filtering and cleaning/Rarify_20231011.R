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

pseq <- Marissa_MU42022_unfiltered

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

#Remove chloroplast, mito, archaea, chimera ----

pseq <- subset_taxa(pseq,Class!="c__Chloroplast")

pseq <- subset_taxa(pseq,Order!="o__Mitochondria")
pseq <- subset_taxa(pseq,Family!="Mitochondria")

pseq <- subset_taxa(pseq,Kingdom!="k__Archaea")
pseq <- subset_taxa(pseq,Kingdom!="Archaea")

pseq <- subset_taxa(pseq,Order!="Chloroplast")

#Check if any chloro, ,mito, or achaeae
tax_levels <- colnames(pseq@tax_table)
chloroplast_archaea_mitochondria <- c("chloroplast", "archaea", "mitochondria")

any_contain <- sapply(chloroplast_archaea_mitochondria, function(term) {
  any(grepl(term, tolower(pseq@tax_table[, tax_levels]), ignore.case = TRUE))
})

any(any_contain)

#If true still contains bad stuff
#if false, move on

rank_names(pseq)

#Remove low prev ----

plot(sort(taxa_sums(x1), TRUE), type="h", ylim=c(0, 10000))

x1 = prune_taxa(taxa_sums(pseq) > 200, pseq) 
x2 = prune_taxa(taxa_sums(pseq) > 500, pseq) 
x3 = prune_taxa(taxa_sums(pseq) > 1000, pseq)

##using x2 mb2021 and x1 for MU42022

library(microbiome)
summarize_phyloseq(pseq)
summarize_phyloseq(x2)


*************************************************
#rarify - removing abundance code with below code does not seem to work...not sure why. Ignore for now.

sequence_abundance <- taxa_sums(pseq)

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

#Check if data needs rarefying
#Rarefication curves ----

#pseq <- x2
pseq <- x1

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
#mu42022: Remove samples T15-S-1, also note, algal samples have very few reads after chloroplasts removed but we will keep them

pseq <- subset_samples(pseq, !Sample.ID %in% c("T15-S-1"))
pseq <- subset_samples(pseq, !Library_Name %in% c("F4L18", "T10r3", "T9r2"))

#saving filtered but not rarefied pseq object for mb2021 project
saveRDS(pseq, "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/mb2021_filtered_NOT_rarefied.rds")
saveRDS(pseq, "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/MU42022_filtered_NOT_rarefied.rds")

otu.rare = otu_table(pseq)
otu.rare = as.data.frame(otu.rare)
sample_names = rownames(otu.rare)

#Generate rarefaction curve, rarefaction curve could be used to determined whether the sequencing depth cover microbial diversity of the sample.
# we will use vegan rarecurve 
library(vegan)

otu.rarecurve <- rarecurve(otu.rare, step = 10000, label = FALSE, xlim = c(0, 100000))

raref.curve <- rarecurve(otu.rare, label = TRUE, ylab = "ASV Count")


#If rarefying:
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

pseq.rel <- microbiome::transform(pseq, "compositional")



pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq_fam <- microbiome::aggregate_rare(ps1, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(ps1, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(ps1, level = "Genus", detection = 50/100, prevalence = 70/100)

##Get a quick look at the data (plot ordination)

set.seed(4235421)

ord <- ordinate(pseq.rel, "MDS", "bray")

#plot MDS/PcoA - can set "colour" and "shape" for any of your variableas
#geompoint controls data point size on plot

plot_ordination(ps1, ord, color = "Family.1", shape = "Percent.Fouled") + geom_point(size = 4) + scale_shape_binned()

plot_ordination(ps1, ord, color = "Factor") + geom_point(size = 4)

plot_ordination(ps1, ord, color = "Family.1") + geom_point(size = 4)

plot_ordination(pseq, ord, color = "Treatment", shape = "Age") + geom_point(size = 4) + facet_wrap(~Age)


p <- plot_ordination(pseq, ord, color = "Factor") + geom_point(size = 4)

