#Sept 27th, 2023
###Generate microbiome heatmaps with top contributing phyla
#source for below: https://rstudio-pubs-static.s3.amazonaws.com/330760_8bba830836324bf6b100d4e76f49e3d2.html#other_visualizations 

#Quick lib ----

library("devtools")
library(phyloseq)
library(microbiome)

#Load data ----


pseq <-  Marissa_MU42022_rare



#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#create objects


pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq, level = "Genus", detection = 50/100, prevalence = 70/100)

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")



###NEXT SECTION - HEATMAPS (TOP ## OTU'S)

#Sort the OTUs by abundance and pick the top 20
top20OTU.names = names(sort(taxa_sums(pseq.fam.rel), TRUE)[1:20])

#Cut down the physeq.tree data to only the top 10 Phyla
top20OTU = prune_taxa(top20OTU.names, pseq.fam.rel)

top20OTU

plot_heatmap(top20OTU)

plot_heatmap(top20OTU, sample.label="Age", sample.order="Age")

plot_heatmap(top20OTU, sample.label="Age", sample.order="Age", taxa.label="Phylum", taxa.order="Family", low="white", high="purple", na.value="grey")

plot_heatmap(top20OTU, "NMDS", "bray", title="Bray-Curtis")


plot_heatmap(top20OTU, sample.label="Tank_treatment", sample.order="Tank_treatment")


plot_heatmap(top20OTU, "NMDS", "bray", title="Bray-Curtis", ample.label="Tank_treatment", sample.order="Tank_treatment")


##I want to look at heatmap of only larval samples

# Exclude specific data from top20OTU object

top20OTU = prune_taxa(top20OTU.names, pseq.fam.rel)

top20OTU

excluded_samples <- c("T10r1", "T10r2", "T10r3", "T11r1", "T11r3", "T12r1", "T12r2", "T12r3", "T13r1","T13r2","T13r3", "T14r1", "T14r2", "T15r1", "T15r2", "T16r1", "T16r2", "T16r3", "T1r1", "T1r2", "T1r3", "T2r1", "T2r3", "T3r1", "T3r2", "T3r3", "T4r1", "T4r2", "T9r1", "T9r3", "T11r2", "T14r3", "T2r2")

top20OTU <- prune_samples(!(sample_names(top20OTU) %in% excluded_samples), top20OTU)

plot_heatmap(top20OTU, sample.label="Tank_treatment", sample.order="Tank_treatment")
