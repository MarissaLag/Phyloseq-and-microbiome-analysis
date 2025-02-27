#Data processing steps from Parker
#Using Cumulative sum scaling (CSS) as normalization step
#Do not have to rarefy data (optional)

#This technique adjusts counts based on the cumulative sum, focusing on ensuring consistent summation across the dataset.

# NORMALIZE WITH CSS

## BiocManager::install("metagenomeSeq")
#install.packages("phyloseq")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("metagenomeSeq")

#remotes::install_github("HCBravoLab/metagenomeSeq")

# install.packages("metagenomeSeq")
# install.packages("zCompositions")
library(metagenomeSeq)
library(phyloseq)
library(zCompositions)

#Data (physeq object)

pseq <- MU42022_filtered_Oct92024

pseq <- PB2023_spat_filtered_not_rarefied

pseq <- PB2023_spat_limited_10X_reads

pseq <- mb2021_filtered_NOT_rarefied

#Preform zComposition (zero imputation) before CSS scaling 

# otu_matrix <- as.matrix(otu_table(pseq))
# 
# otu_imputed <- cmultRepl(otu_matrix, 
#                          method = "CZM")
# 
# otu_table_imputed <- otu_table(otu_imputed, taxa_are_rows = FALSE)
# 
# pseq_imputed <- phyloseq(otu_table_imputed, sample_data(pseq), tax_table(pseq))

#If samples removed, must remove from sample_data

# Get sample names from both phyloseq objects
# original_samples <- sample_names(pseq)
# imputed_samples <- sample_names(pseq_imputed)
# 
# # Find samples that were removed
# removed_samples <- setdiff(original_samples, imputed_samples)
# print(removed_samples)
# 
# # Subset sample_data to keep only samples in pseq_imputed
# sample_data(pseq_imputed) <- sample_data(pseq_imputed)[imputed_samples, ]
# 
# pseq <- pseq_imputed

# Convert from filtered phyloseq object to metagenomeseq object:
metaSeqObject1=phyloseq_to_metagenomeSeq(pseq)
p=cumNormStatFast(metaSeqObject1)
metaSeqObject_CSS=cumNorm(metaSeqObject1, p) # CSS normalization function
metaSeqObject_CSS_filt=filterData(metaSeqObject_CSS, depth=1000) #Excludes low reads, minimum depth is 1000.
# With 1000 minimum depth, 5 samples were removed
seq.asv.css=data.frame(MRcounts(metaSeqObject_CSS_filt, norm=TRUE, log=TRUE)) # Convert back to an asv table, now css corrected.

# Fixes column names in the ASV table to match our metadata by erasing the X character at the beginning. (Something like this may or may not be needed for your data, you just need the sample names to match before remaking the phyloseq object)
names(seq.asv.css) = gsub(pattern = "...", replacement = "---", x = names(seq.asv.css))
# Replace "..." with "---" in the column names of the data frame
colnames(seq.asv.css) <- gsub("\\.\\.\\.", "---", colnames(seq.asv.css))

seq.asv.css <- t(seq.asv.css)

# Remove samples that were filtered out from the metadata file before remaking phyloseq object.
# sample.names2=colnames(seq.asv.css) #Pull the sample names from the sequence table.
# deleted=print(setdiff(sample.names,sample.names2)) #Creates an object listing the sample names that were removed during normalization.
# metadata=metadata[!row.names(metadata) %in% deleted,] #Delete samples that were removed during normalization from metadata

# Recreate phyloseq object.
taxa = pseq@tax_table
metadata = pseq@sam_data

physeq.sub.arch=phyloseq(otu_table(seq.asv.css, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa))
physeq.sub.arch
# If sample names are column names in the sequence table, taxa_are_rows=TRUE
# If samples names are row names, taxa_are_rows=FALSE

saveRDS(physeq.sub.arch, file = "mb2021_filtered_NOT_rarefied_CSS.rds")


