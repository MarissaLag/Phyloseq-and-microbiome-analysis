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

install.packages("metagenomeSeq")
library(metagenomeSeq)
library(phyloseq)

#Data (physeq object)

pseq <- MU42022_filtered_Oct92024

# Convert from filtered phyloseq object to metagenomeseq object:
metaSeqObject1=phyloseq_to_metagenomeSeq(pseq)
p=cumNormStatFast(metaSeqObject1)
metaSeqObject_CSS=cumNorm(metaSeqObject1, p) # CSS normalization function
metaSeqObject_CSS_filt=filterData(metaSeqObject_CSS, depth=1000) #Excludes low reads, minimum depth is 1000.
# With 1000 minimum depth, 5 samples were removed
seq.asv.css=data.frame(MRcounts(metaSeqObject_CSS_filt, norm=TRUE, log=TRUE)) # Convert back to an asv table, now css corrected.

# Fixes column names in the ASV table to match our metadata by erasing the X character at the beginning. (Something like this may or may not be needed for your data, you just need the sample names to match before remaking the phyloseq object)
names(seq.asv.css) = gsub(pattern = "X", replacement = "", x = names(seq.asv.css))

# Remove samples that were filtered out from the metadata file before remaking phyloseq object.
# sample.names2=colnames(seq.asv.css) #Pull the sample names from the sequence table.
# deleted=print(setdiff(sample.names,sample.names2)) #Creates an object listing the sample names that were removed during normalization.
# metadata=metadata[!row.names(metadata) %in% deleted,] #Delete samples that were removed during normalization from metadata

# Recreate phyloseq object.
physeq.sub.arch=phyloseq(otu_table(seq.asv.css, taxa_are_rows=TRUE), sample_data(metadata), tax_table(taxa))
physeq.sub.arch
# If sample names are column names in the sequence table, taxa_are_rows=TRUE
# If samples names are row names, taxa_are_rows=FALSE
# I want to convert back to a phyloseq object so I can use the microbiome package to run an NMDS