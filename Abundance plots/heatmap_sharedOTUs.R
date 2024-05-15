#Heat map of top ASVs
#trying to replicate heatmap done in Domin et al., 2023
#source: https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_extract_otus.R

devtools::install_github("vmikk/metagMisc")
install.packages("metagMisc")

library("devtools")
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)
library(metagMisc)
library(pheatmap)

#Select ASVs "which were constantly detected in at least one of the four 
#different timepoints of recolonization across all samples (minimum relative abundance was 0.005%)" (Domin et al., 2023)

#For James' experiment, select core ASVs within treatments Control and Antibiotics

pseq <-`Filtered_Rarified_MU42022_23-12-13`

#Exclude factors

pseq_filt <- subset_samples(pseq, Age %in% c("Day 01", "Day 03", "Day 06", "Day 15"))

pseq_filt <- subset_samples(pseq_filt, Treatment %in% c("Control", "Antibiotics"))

OTU = pseq_filt@otu_table
Tax = pseq_filt@tax_table
Metadata = pseq_filt@sam_data
Tree = pseq_filt@phy_tree

View(OTU)
samplenames = Metadata$Treatment



View(Metadata)
#Taxa shared between all samples samples

ps <- phyloseq_extract_shared_otus(pseq_filt, samp_names = pseq_filt)
ps <- phyloseq_extract_shared_otus(pseq_filt, samp_names = c("T1-1", "T1-15", "T1-3", "T1-6", "T13-1", "T13-15", "T13-3", "T13-6", "T16-1-r1", "T16-1-r2", "T16-15-r1", "T16-15-r2", "T16-15-r3", "T16-3-r1", "T16-3-r2", "T16-6-r1", "T16-6-r2", "T19-1", "T19-3", "T19-6", "T7-1", "T7-15", "T7-3", "T7-6" ))

pseq_filt@otu_table[1:5,1:5]
head(pseq_filt@sam_data)
pseq_filt@sam_data[pseq_filt@sam_data$Treatment=="Control", "Sample.ID"]

setdiff()
intersect()

pseq_filt

# Make a subset of pseq_filt with only Control samples, get all taxa IDs

# 
?subset_samples()

###Andy Try
#source: https://microbiome.github.io/tutorials/Core.html
ps10<-merge_samples(pseq_filt, "Treatment", fun= mean)
ps10@otu_table[,1:10]
?core_members
Core_trial<-core(ps10, detection=0.5/100, prevalence = 90/100)
View(Core_trial@otu_table)
Core_trial

Core_trial <- psmelt(Core_trial)
str(Core_trial)

#create list of ASV names

otu_names_list <- Core_trial$OTU
otu_names_list




#recall pseq_filt object
pseq_filt

#prune taxa to only include taxa on list (shared between treatments)

otu_names_list <- Core_trial$OTU

prune_pseq_filt<- prune_taxa(otu_names_list, pseq_filt)

prune_pseq_filt



#713 taxa remaining, lets reduce it to 100 taxa to make heat map clearer
#top 300 most abundant taxa
prune_pseq_filt <- prune_taxa(names(sort(taxa_sums(prune_pseq_filt),TRUE)[1:100]), prune_pseq_filt)
prune_pseq_filt@sam_data

plot_heatmap(prune_pseq_filt, sample.label="Age")

sample_data(prune_pseq_filt) <- sample_data(prune_pseq_filt)[, !(colnames(sample_data(prune_pseq_filt)) %in% "new_column_name")]

Treat_Age <- c("C - 1dpf", "C - 15dpf", "C - 3dpf", "C - 6dpf", "C - 1dpf", "C - 15dpf", "C - 3dpf", "C - 6dpf", "A - 1dpf", "A - 1dpf", "A - 15dpf", "A - 15dpf", "A - 15dpf", "A - 3dpf", "A - 3dpf", "A - 6dpf", "A - 6dpf", "C - 1dpf", "C - 3dpf", "C - 6dpf", "C - 1dpf", "C - 15dpf", "C - 3dpf", "C - 6dpf")
prune_pseq_filt@sam_data$Treat_Age <- Treat_Age

prune_pseq_filt_comb <-merge_samples(prune_pseq_filt, "Treat_Age", fun= mean)
View(prune_pseq_filt_comb@sam_data)

#plot treatments seperately
prune_pseq_filt <- subset_samples(prune_pseq_filt, Treatment %in% c("Control"))
sample_order <- c("Day 01", "Day 03", "Day 06", "Day 15")
plot_heatmap(prune_pseq_filt, sample.label="Age", sample.order = "Age", method = "NMDS", distance = "bray")
p <- plot_heatmap(prune_pseq_filt, sample.label="Age", sample.order = "Age", method = "NMDS", distance = "bray", facet_grid("Treatment"))
p <- plot_heatmap(prune_pseq_filt_comb, sample.label="Age", sample.order = "Age", method = "NMDS", distance = "bray", facet_grid("Treatment"))
p + facet_wrap("Treatment")


order <- c("C - 1dpf", "A - 1dpf", "C - 3dpf", "A - 3dpf", "C - 6dpf", "A - 6dpf", "C - 15dpf", "A - 15dpf")
order <- c("C - 1dpf", "C - 3dpf", "C - 6dpf", "C - 15dpf", "A - 1dpf", "A - 3dpf", "A - 6dpf", "A - 15dpf")

plot_heatmap(prune_pseq_filt_comb, sample.order = order)


#Only map ASVs found to be signif different between Age groups

asv_list_all <- dimnames(subset_data)[[1]]

age_ASVs <- prune_taxa(asv_list, pseq_filt)


#convert to compositional data

age_ASVs.rel <- microbiome::transform(age_ASVs, "compositional")

psmelt_age_ASVs <- psmelt(age_ASVs.rel)

View(psmelt_age_ASVs)

plot_heatmap(age_ASVs, sample.order = order)

plot_heatmap(age_ASVs.rel, sample.label = "Sample.ID",
           high="#66CCFF", low="#000033")

plot_heatmap(age_ASVs, sample.label = "Age")

#pheatmaps ----
install.packages("dendextend")
library(dendextend)

#create dendrogram of variable (e.g., genes or ASVs)

my_hclust_ASV <- hclust(dist(psmelt_age_ASVs), method = "complete")

as.dendrogram(my_hclust_ASV) %>%
  plot(horiz = TRUE)

#cut tree into clusters - here splitting into 2 clusters

my_ASV_col <- cutree(tree = as.dendrogram(my_hclust_ASV), k = 4)

#rename to cluster name

#if 2 clusters (i.e., k = 2) selected
my_ASV_col <- data.frame(cluster = ifelse(test = my_ASV_col == 1, yes = "cluster 1", no = "cluster 2"))

#if more clusters used
my_ASV_col <- data.frame(cluster = factor(my_ASV_col, labels = c("cluster 1", "cluster 2", "cluster 3", "cluster 4")))

head(my_ASV_col)


#add some column annotations and create the heatmap.
row.names() <- colnames(inv_F$str)
pheatmap(psmelt_age_ASVs, annotation_row = my_ASV_col, annotation_col = my_sample_col)

pheatmap(subset_data, 
         annotation_row = my_ASV_col, 
         show_rownames = FALSE,
         angle_col = 0,
         cluster_cols = FALSE)
