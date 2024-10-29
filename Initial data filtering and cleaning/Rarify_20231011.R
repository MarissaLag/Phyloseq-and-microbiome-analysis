##testing MDJ.rds data

#Load packages ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

library("devtools")
library(phyloseq)
library(microbiome)
library(phyloseq)


#Load un-rarefied data ----

pseq<- readRDS("Marissa_.rds")

pseq <- Marissa_mb2021_filtered_20240203

pseq <- Marissa_MU42022_unfiltered

pseq <- Marissa_mb2021_unfiltered

pseq


#SMK project - seperate into different RDS objects
# Subset the data for each project
SMK_Sam <- subset_samples(SMK_2024, Project == "Sam")
SMK_Marissa <- subset_samples(SMK_2024, Project == "PB2023")
SMK_Korrina <- subset_samples(SMK_2024, Project == "Korrina")

# Save each subset as an .rds file
saveRDS(SMK_Sam, "SMK_Sam.rds")
saveRDS(SMK_Marissa, "PB2023.rds")
saveRDS(SMK_Korrina, "SMK_Korrina.rds")

pseq <- SMK_Marissa

#removing samples

pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")
pseq <- subset_samples(pseq, Age %in% c("Spat"))

pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))

pseq <- subset_samples(pseq, Age %in% c("Spat"))

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

# Create table, number of features for each phyla
#Checking for NA phyla - these are often sequence artifacts and should be removed
table(tax_table(pseq)[, "Phylum"], exclude = NULL)

#If need to remove: 
#pseq <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove low prev ----
#remove less than 0.5% total abundance - 5234 taxa * 0.005 = ~30 reads minimum (0.5%)
plot(sort(taxa_sums(x2), TRUE), type="h", ylim=c(0, 10000))
# 
# x0 = prune_taxa(taxa_sums(pseq) > 30, pseq) 
x1 = prune_taxa(taxa_sums(pseq) > 200, pseq) 
x2 = prune_taxa(taxa_sums(pseq) > 500, pseq) #use for PB2023
x3 = prune_taxa(taxa_sums(pseq) > 1000, pseq)

summarize_phyloseq(pseq)
summarize_phyloseq(x2)

pseq <- x2

##using x2 mb2021 and x1 for MU42022

#Removing low prevalent taxa

#rarify - removing abundance code with below code does not seem to work...not sure why. Ignore for now.

#Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(pseq),
               MARGIN = ifelse(taxa_are_rows(pseq), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(pseq),
                    tax_table(pseq))

#Check for low prevelance phyla and remove
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#mb2021: Armatimonsdota and Sumerlaetoa only in one sample so will remove
#mb2021: Entotheonella, Gemmatimondota, Patescibacteria only in 4, 5, and 3 samples so will remove
# Define phyla to filter
filterPhyla = c("Acidobacteriota", "Chloroflexi")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(pseq, !Phylum %in% filterPhyla)
ps1


# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(pseq),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
  
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(pseq)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps1)

ps2
plot(sort(taxa_sums(ps2), TRUE), type="h", ylim=c(0, 10000)) #still has tail

#Agglomerate taxonomic redundancies
# How many genera would be present after filtering?
length(get_taxa_unique(ps3, taxonomic.rank = "Genus"))
## [1] 363
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
#362 - since hardly different, will skip this



#Abundance transformations
#It is usually necessary to transform microbiome count data to account for differences in library size, variance, scale, etc
plot_abundance = function(pseq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(pseq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Treatment",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
#The transformation in this case converts the counts from each sample into their frequencies, often referred to as proportions or relative abundances
# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(ps2, function(x){x / sum(x)})
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)

#Look for taxa will bimodal abundance
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)

#Depending on data, may not always have to rarefy read counts
#if rarefication curves look good, don't rarefy. If there are one or two outliers, remove them and don't rarefy
#Also Jadi normalizes her data as well (even if rarefied or not) with DESEq2 negative binomial dist'd model
#Check if data needs rarefying
#Rarefication curves ----

pseq <- x2 #899 taxa for mb2021
pseq <- x1
pseq <- x0
pseq <- ps2 #2609 taxa (filter >5%)
#check reads

library(ggplot2) 
library(data.table)

# Check read count
readcount = data.table(as(sample_data(pseq), "data.frame"),
                         TotalReads = sample_sums(pseq), 
                         keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID")

ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")

View(readcount[order(readcount$TotalReads), c("SampleID", "TotalReads")])

readcount[order(readcount$TotalReads), c("SampleID", "TotalReads")]

#mb2021: Remove samples F4L18, T10r3, T9r2 (did not work at all)
#mu42022: Remove samples T15-S-1, also note, algal samples have very few reads after chloroplasts removed but we will keep them
#PB2023: Remove T7-2-S (<200 reads)

pseq <- subset_samples(pseq, !Sample.ID %in% c("T7-2-S"))
pseq <- subset_samples(pseq, !Library_Name %in% c("F4L18", "T10r3", "T9r2"))

#saving filtered but not rarefied pseq object for mb2021 project
saveRDS(pseq, "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/mb2021_filtered_NOT_rarefied.rds")
saveRDS(pseq, "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/MU42022_filtered_Oct92024.rds")
saveRDS(pseq, "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/mb2021_filteredwSpat_only_rarefied_June2024.rds")
saveRDS(pseq, "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/MU42022_spat_unrarefied.rds")


otu.rare = otu_table(pseq)
otu.rare = as.data.frame(otu.rare)
sample_names = rownames(otu.rare)

#Generate rarefaction curve, rarefaction curve could be used to determined whether the sequencing depth cover microbial diversity of the sample.
# we will use vegan rarecurve 
library(vegan)

otu.rarecurve <- rarecurve(otu.rare, step = 10000, label = TRUE, xlim = c(0, 10000))

raref.curve <- rarecurve(otu.rare, label = TRUE, ylab = "ASV Count")

#If not rarefying - DESeq2 transformation ----
#counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") #download deseq

BiocManager::install("DESeq2") #Install deseq

BiocManager::install("DESeq2", force = TRUE)

BiocManager::install("GenomeInfoDb")
BiocManager::install("DESeq2")

library(DESeq2)
library(phyloseq)

DeSeq <- phyloseq_to_deseq2(pseq, ~ Treatment) #convert phyloseq to deseq object

DeSeq2 <- DESeq(DeSeq)

# # countData = ASV counts with samples in columns  
# # colData = Sample meta data (factors, etc)  
# # design = the right hand side of a GLM formula

ASV_deseq <- DESeqDataSetFromMatrix(countData = t(),
                                    colData = fact2,
                                    design = ~ Treatment)

# positive counts normalisaton (accounts for zero-inflation)
ASV_deseq <- estimateSizeFactors(DeSeq, type = "poscounts")
sizeFactors(ASV_deseq)

# size-factor corrected data are calculated by dividing the raw counts by the sample size factor and adding 0.5 to correct for the zeros
ASV_deseq_norm <- sapply(row.names(ASV_deseq), function(x){
  plotCounts(ASV_deseq, x, "Day", returnData = TRUE, normalized = TRUE)$count
  
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
write.csv(ASV_deseq_norm, "~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files//normalised_ASV_table.csv")

#Make new pseq object with normalized ASV table
mb2021_filtered_NOT_rarefied <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/mb2021_filtered_NOT_rarefied.rds")

#Check format
View(mb2021_filtered_NOT_rarefied@otu_table)
View(normalised_ASV_table)

#correct to make forst column as row names

normalised_ASV_table <- as.data.frame(normalised_ASV_table)

rownames(normalised_ASV_table) <- normalised_ASV_table[, 1]

# Remove the first column from the data frame
normalised_ASV_table <- normalised_ASV_table[, -1]

# Create a new otu_table object
new_otu_table <- otu_table(normalised_ASV_table, taxa_are_rows = FALSE)

# Replace the old OTU table with the new one
otu_table(mb2021_filtered_NOT_rarefied) <- new_otu_table

# Verify the replacement
View(mb2021_filtered_NOT_rarefied@otu_table)
#Save as different RDS file
saveRDS(mb2021_filtered_NOT_rarefied, file = "mb2021_filtered_NOT_rarefied_normalized.rds")


#If rarefying:
##rarify data to make sequences an even read depth - selecting read depth of 10,000 = any samples with fewer than 10,000 total reads will be removed, all samples will be equalized to 5000 reads for Denman's samples, 10,000 for Marissa/James'

Rare <-rarefy_even_depth(pseq, sample.size= 5128)

#Marissa/James = 17 samples removed because they contained fewer reads than `sample.size' - first 5 reads are T13-1,T13-2-2,T14-2-2,T15-S-1,T16-10-r3

##check if seq depth is 10,000 or 5,000

#mb2021 rarefied to 3000

sample_depths <- sample_sums(Rare)

print(sample_depths)

##yes, all samples have 10,000 or 5,000 reads now.

pseq <- Rare

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

plot_ordination(pseq, ord, color = "Treatment", shape = "Family") + geom_point(size = 4) + scale_shape_binned()

