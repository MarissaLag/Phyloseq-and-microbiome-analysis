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
miseq_path <- "R1:R2fastq_files/"
list.files(miseq_path)

# Sort ensures forward/reverse reads are in same order. notice the pattern (two different reads, Forward and Reverse)
fnFs <- sort(list.files(miseq_path, pattern="_R1.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2.fastq"))
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
plotQualityProfile(fnFs[1:8])

plotQualityProfile(fnRs[1:8])

#for mb2021 - forward reads good, reverse good (above QS of 30) until 200bp - may want to remove reverse
####Quality for Reverse is bad. So we are only doing the Forward 

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

#### We can use this data top say where within the sequence to trim the data. "Here, the forward reads maintain high quality throughout, while the quality of the reverse reads drops significantly at about position 160. Therefore, we choose to truncate the forward reads at position 200, and the reverse reads at position 160. We also choose to trim the first 25 nucleotides of each read based on empirical observations across many Illumina datasets that these base positions are particularly likely to contain pathological errors."
length(fnFs)
length(fnRs)
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150, 200),
                      maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft = 10,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)



#Make ASVs ----
derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
#names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)

#errR <- learnErrors(filtRs, multithread=TRUE)
 
#plotErrors(errR)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#Check ASVs - if many unknowns, you may want to change your trimming parameters (make them larger) to provide more info

##Still have many unique sequences - (about 1/3rd of reads unique seqs) - change to 100 or 120bp region? -> actually, after looking at taxonomy_alldata.csv, looks okay

#If using both forward and reverse
#mergers <- mergePairs(dadaFs, derepFs) 

seqtabAll <- makeSequenceTable(dadaFs[!grepl("Mock", names(dadaFs))])
table(nchar(getSequences(seqtabAll)))

seqtabNoC <- removeBimeraDenovo(seqtabAll)

fastaRef <- "./silva_nr99_v138.1_train_set.fa"

taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))

write.csv(taxTab, "taxonomy_alldata_2.csv")


#ASVs to Taxa info ----
# now replace the long ASV names (the actual sequences) with human-readable names
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtabNoC)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names

write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.txt", sep="\t", quote=F, row.names=F, col.names=F)

library(seqinr)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "16s_ASV_sequences_all.fasta") #save sequences with new names in fasta format


#Phylogenetic tree ----

###Try 
Object1<- cbind(ASV.num, taxTab)

#IMPORTANT: sanity checks
colnames(seqtabNoC) == ASV.seq #only proceed if this tests as true for all elements
row.names(taxTab) == ASV.seq #only proceed if this tests as true for all elements

#rename your ASVs in the taxonomy table and sequence table objects
colnames(seqtabNoC) <- ASV.num
row.names(taxTab) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(seqtabNoC),seqtabNoC),"sequence_table.16s.all_merged.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(taxTab),taxTab),"taxonomy_table.16s_all_merged.txt", row.names=FALSE, quote=F, sep="\t")


#### Phylogenetic tree 
library(phangorn)
Fast1<- readDNAStringSet(file= "./16s_ASV_sequences.fasta",format = "fasta")

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
 
meta<-import_qiime_sample_data("Metadata_Marissa_only.txt")

ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), tax_table(taxTab))

ps1 <-merge_phyloseq(ps,meta)
ps1

##MDJ.rds has James and Marissa's samples only
#Marissa_MU42022 has marissa's samples only
#Denman has Denman's only

saveRDS(ps1, file= "Marissa_MU42022.rds")

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

pseq <- Denman_samples


##### Filter out low frequency and low abundance reads. 
pseq_filter<-phyloseq_filter_prevalence(pseq, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR", abund.type = "total")


phyloseq_filter_prevalence()

######Rare_faction 

Rare_10000<-rarefy_even_depth(pseq_filter, sample.size= 10000)







