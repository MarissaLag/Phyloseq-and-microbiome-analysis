#Assessing normality and zeros in microbiome data
#Adding visuals to help assess

#Load packages ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

install.packages("phyloseq")
install.packages("devtools")

library("devtools")
library(phyloseq)
library(microbiome)
library(gridExtra)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") #download deseq

BiocManager::install("DESeq2") #Install deseq

BiocManager::install("DESeq2", force = TRUE)

BiocManager::install("GenomeInfoDb")
BiocManager::install("DESeq2")
BiocManager::install("IRanges", ask = FALSE, force = TRUE)

library(DESeq2)

pseq <- MU42022_filtered_Oct92024
pseq <- MU42022_filtered_NOT_rarefied

pseq


# Extract OTU table from phyloseq object and ensure it's in matrix format
otu_mat <- as(otu_table(pseq), "matrix")

# Check if OTUs are rows; if not, transpose the matrix
if (!taxa_are_rows(pseq)) {
  otu_mat <- t(otu_mat)
}

# Calculate row sums (total abundance per sample)
sample_sums <- rowSums(otu_mat)

# Plot histogram of total abundance per sample
hist(sample_sums, main = "Total Abundance per Sample", xlab = "Abundance")

#First step is to use a normalization method
#Using Deseq2 - expects count data, so cannot be in relative abundace
#Deseq considered a "scaling method" - although it use for 
#Differential abundance is questioned, it use for zero-inflated data is pretty good

#Cannot have NAs in factor of interest
# Access the sample_data from your phyloseq object
sam_data <- sample_data(pseq)

# Replace NA values in the Treatment column with "Algae"
sam_data$Treatment[is.na(sam_data$Treatment)] <- "Algae"

# Update the sample_data in your phyloseq object
sample_data(pseq) <- sam_data

# Check if the changes were applied
table(sample_data(pseq)$Treatment)

sample_data(pseq)$Treatment <- as.factor(sample_data(pseq)$Treatment)
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))
pseq <- subset_samples(pseq, !Treatment %in% c("Algae"))

#Sanity check
#check if any OTUs are not present in any samples (want false)
any(taxa_sums(pseq) == 0)

#if true
pseq_filtered <- prune_taxa(taxa_sums(pseq) > 0, pseq)
any(taxa_sums(pseq_filtered) == 0)

pseq <- pseq_filtered

# Convert phyloseq to DESeq2 object
#Choose afctor you are largely interested - tells Deseq what factor you are looking at
#So it can account for this while normalizing


DeSeq <- phyloseq_to_deseq2(pseq, ~ Treatment) #convert phyloseq to deseq object
summary(assay(DeSeq))

DeSeq2 <- DESeq(DeSeq)

# Get normalized counts
normalized_counts <- counts(DeSeq2, normalized = TRUE)
# Summary of normalized counts
summary(normalized_counts)

# Get results for Treatment
results_deSeq <- results(DeSeq2)

# View results
head(results_deSeq)

# Optionally filter results for significant taxa
sig_results <- results_deSeq[which(results_deSeq$padj < 0.05), ]
cat("Number of significant taxa:", nrow(sig_results), "\n")

# MA plot
plotMA(results_deSeq, ylim = c(-5, 5), main = "MA Plot")

#Size factors
ASV_deseq <- estimateSizeFactors(DeSeq2, type = "poscounts")

#Each sample's size factor reflects the relative scaling needed to normalize its counts. 
#A size factor of 1 indicates that the sample does not require scaling; 
#values greater than 1 indicate that the counts for that sample are scaled up, while values less than 1 indicate that the counts are scaled down.
sizeFactors(ASV_deseq)

# size-factor corrected data are calculated by dividing the raw counts by the sample size factor and adding 0.5 to correct for the zeros
ASV_deseq_norm <- sapply(row.names(ASV_deseq), function(x){
  plotCounts(ASV_deseq, x, "Treatment", returnData = TRUE, normalized = TRUE)$count
  
})
rownames(ASV_deseq_norm) <- colnames(ASV_deseq)

ASV_deseq_norm

# remove the addition of 0.5 to all entries
ASV_deseq_norm <- ASV_deseq_norm - 0.5
ASV_deseq_norm[1:10, 1:10]

# make dataset with integers
ASV_deseq_norm <- round(ASV_deseq_norm)
ASV_deseq_norm[1:10, 1:10]

# check for zeros across all samples (produced by normalisation)
dim(ASV_deseq_norm[, colSums(ASV_deseq_norm) == 0]) # 3

# check singletons
dim(ASV_deseq_norm[, colSums(ASV_deseq_norm) == 1]) # 5 singletons - keep these


#Compare normalized and non-normalized otu table
ASV_deseq_norm #normalized
OTU <- pseq@otu_table

# Check if row names (OTUs) and column names (samples) match
all(rownames(OTU) == rownames(ASV_deseq_norm)) # TRUE if they match
all(colnames(OTU) == colnames(ASV_deseq_norm)) # TRUE if they match

# You can calculate differences or correlations
comparison <- OTU - ASV_deseq_norm

# View the comparison
head(comparison)


#Put back into pseq object
new_otu_table <- otu_table(ASV_deseq_norm, taxa_are_rows = FALSE)

# Replace the old OTU table with the new one
otu_table(pseq) <- new_otu_table

pseq@otu_table

#Convert abundace into relative abundance/propotional data
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
ps2ra = transform_sample_counts(pseq, function(x){x / sum(x)})
plotBefore = plot_abundance(pseq, "")
plotAfter = plot_abundance(ps2ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)


#Save normalized, relative abundance data 
saveRDS(pseq, "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/MU42022_normalized_relative.rds")

pseq


