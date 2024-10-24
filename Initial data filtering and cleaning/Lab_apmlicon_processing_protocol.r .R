# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html

#Packages ----
### Downloading all the packages to process data in R
setwd("/Users/greent/Desktop/Amplicon_Seq_data_Aug_2023/Marissa_Seqs")

# Clear workspace 
# rm(list=ls())

install.packages("devtools")
library("devtools")
# devtools::install_github("benjjneb/dada2"
#                          #, ref="v1.16"
#                          ) # change the ref argument to get other versions

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")
BiocManager::install("phangorn")
BiocManager::install("dada2")
BiocManager::install("BiocStyle")

install.packages("seqinr")
library("dada2")
library("seqinr")
library("knitr")
library("BiocStyle")


.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")


# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

set.seed(100)

#Load seqs ----
####Tell R where the data is...
miseq_path <- "R1R2_SMK_2024/"
list.files(miseq_path)

# Sort ensures forward/reverse reads are in same order. notice the pattern (two different reads, Forward and Reverse)
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
#fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
#fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)

#Did it work? 
fnFs[1:3]
fnRs[1:3]

#Seq quality ----
#Quality of reads: Most Illumina sequencing data shows a trend of decreasing average quality towards the end of sequencing reads.
#This only shows you the first 2. Which direction is this for? 
plotQualityProfile(fnFs[1:10])

plotQualityProfile(fnRs[50:59])

#for mb2021 - forward reads good, reverse good (above QS of 30) until 200bp - may want to remove reverse
####Quality for Reverse is bad. So we are only doing the Forward - MU42022

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

#### We can use this data top say where within the sequence to trim the data. "Here, the forward reads maintain high quality throughout, 
#while the quality of the reverse reads drops significantly at about position 160. 
#Therefore, we choose to truncate the forward reads at position 200, and the reverse reads at position 160. 
#We also choose to trim the first 25 nucleotides of each read based on empirical observations across many Illumina datasets 
#that these base positions are particularly likely to contain pathological errors."

length(fnFs)
length(fnRs)
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250, 180),
                      maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft = 10,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Plot quality again

plotQualityProfile(filtFs[148:149])

plotQualityProfile(filtRs[70:80])

#Make ASVs ----
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#If getting error "files do not exist" from above, run this code:
file_check <- file.exists(filtFs) 
false_indices <- which(!file_check)
missing_files <- filtFs[false_indices]
print(missing_files)
#SMK project- samples: F3-NPC-PCR and F3-PC-PCR had no reads so were removed

filtFs_cleaned <- filtFs %>%
  setdiff(c("R1R2_SMK_2024//filtered/F3-NPC-PCR_F_filt.fastq.gz", 
            "R1R2_SMK_2024//filtered/F4-PC-PCR_F_filt.fastq.gz"))

filtRs_cleaned <- filtRs %>%
  setdiff(c("R1R2_SMK_2024//filtered/F3-NPC-PCR_R_filt.fastq.gz", 
            "R1R2_SMK_2024//filtered/F4-PC-PCR_R_filt.fastq.gz"))

#rename to original file and run code from above:
filtFs <- filtFs_cleaned
filtRs <- filtRs_cleaned

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

#If removed any files (no reads) will have to update sampleNames object so lengths match
sampleNames_new <- sampleNames %>%
  setdiff(c("F3-NPC-PCR", 
            "F4-PC-PCR"))

sampleNames <- sampleNames_new

#learn error rates - long step
errF <- learnErrors(filtFs, multithread=TRUE)

#Got this message: Warning message:
#In browseURL(paste0("http://127.0.0.1:", port, "/library/", pkgname,  :
#closing unused connection 3 (R1R2_SMK_2024//filtered/120_F_filt.fastq.gz)

errR <- learnErrors(filtRs, multithread=TRUE)
 
plotErrors(errF)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#Check ASVs - if many unknowns, you may want to change your trimming parameters (make them larger) to provide more info

##Still have many unique sequences - (about 1/3rd of reads unique seqs) - change to 100 or 120bp region? -> actually, after looking at taxonomy_alldata.csv, looks okay

#If using both forward and reverse
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs) 

seqtabAll <- makeSequenceTable(dadaFs[!grepl("Mock", names(dadaFs))])
table(nchar(getSequences(seqtabAll)))

seqtabNoC <- removeBimeraDenovo(seqtabAll)

#SMK2024 project - using newest SILVA database (138.2 versus 138.1)
fastaRef <- "./silva_nr99_v138.1_train_set.fa"

taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))

write.csv(taxTab, "taxonomy_alldata_2_SMK.csv")


#ASVs to Taxa info ----
# now replace the long ASV names (the actual sequences) with human-readable names
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtabNoC)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names

write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping_SMK.txt", sep="\t", quote=F, row.names=F, col.names=F)

library(seqinr)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "16s_ASV_sequences_all_SMK.fasta") #save sequences with new names in fasta format



#Phylogenetic tree ----

###Try 
Object1<- cbind(ASV.num, taxTab)

#IMPORTANT: sanity checks
colnames(seqtabNoC) == ASV.seq #only proceed if this tests as true for all elements -true
row.names(taxTab) == ASV.seq #only proceed if this tests as true for all elements -true

#rename your ASVs in the taxonomy table and sequence table objects
colnames(seqtabNoC) <- ASV.num
row.names(taxTab) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(seqtabNoC),seqtabNoC),"sequence_table.16s.all_merged_SMK.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(taxTab),taxTab),"taxonomy_table.16s_all_merged_SMK.txt", row.names=FALSE, quote=F, sep="\t")


#### Phylogenetic tree 
library(phangorn)

install.packages("BiocManager")
BiocManager::install("Biostrings")
detach("package:Rsamtools", unload = TRUE)
detach("package:GenomicRanges", unload = TRUE)
# Repeat for other packages as necessary
BiocManager::install("S4Vectors")

library(Biostrings)

Fast1<- readDNAStringSet(file= "./16s_ASV_sequences_SMK.fasta",format = "fasta")

seqs <- getSequences(Fast1)
names(Fast1) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)


plot(fitGTR)

#RDS phyloseq object ----
 
meta<-import_qiime_sample_data("Meta_all.txt")

ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), tax_table(taxTab))

ps1 <-merge_phyloseq(ps,meta)
ps1

#Remove Na columns in dataframe
sample_data_df <- data.frame(sample_data(ps1))
sample_data_df <- sample_data_df[, !grepl("^X", colnames(sample_data_df))]
sample_data(ps1) <- sample_data(sample_data_df)

##MDJ.rds has James and Marissa's samples only
#Marissa_MU42022 has marissa's samples only
#Denman has Denman's only

saveRDS(ps1, file= "SMK_2024.rds")

MDS<- ordinate(ps1, method = "NMDS", distance = "bray", weighted = TRUE)
MDS_Bray<- plot_ordination(ps1, MDS, color = "Treatment")

MDS1<- ordinate(ps1, method = "MDS", distance = "sor")
MDS_Bray1<- plot_ordination(ps1, MDS1, color = "Treatment")+ theme_bw()


Uni_w<- ordinate(ps1, method = "MDS", distance = "unifrac", weighted = TRUE)
Ph_w<- plot_ordination(ps1, Uni_w, color = "Treatment") + theme_bw() + title("Unifrac weighted")

Uni_u<- ordinate(ps1, method = "MDS", distance = "unifrac", weighted = FALSE)
Ph_u<- plot_ordination(ps1, Uni_u, color = "Treatment")
grid.arrange(Ph_w,Ph_u,MDS_Bray, MDS_Bray1)


#Rarefying data ----

#Depending on data, may not always have to rarefy
#if rarefication curves look good, don't rarefy. If there are one or two outliers, remove them and don't rarefy
#Also Jadi normalizes her data as well (even if rarefied or not) with DESEq2 negative binomial dist'd model
#see Rarify_20231011.R for details on how to check this. 

#IF rarefying: Remove low freq taxa and standardize read count
##start pre-processing

install.packages("remotes")
remotes::install_github("vmikk/metagMisc")

#Libraries to load:
library(phyloseq)
library(microbiome)
install.packages("data.table")
library(data.table)
library(vegan)
library("metagMisc")
library(metagMisc)

#To standardize our sampling we rarefy the data. Essentially we make an arbitrary cutoff of how many samples we will examine.  We use the rarefaction, plus the list of sequences per sample to make this call (check out object "sdt" that we made earlier. In this dataset we will do 10,000. 




##### Filter out low frequency and low abundance reads. 
pseq_filter<-phyloseq_filter_prevalence(pseq, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR", abund.type = "total")


phyloseq_filter_prevalence()

######Rare_faction 

Rare_10000<-rarefy_even_depth(pseq_filter, sample.size= 10000)







