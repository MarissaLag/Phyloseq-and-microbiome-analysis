#Sept 27th, 2023
###Generate microbiome heatmaps with top contributing phyla
#source for below: https://rstudio-pubs-static.s3.amazonaws.com/330760_8bba830836324bf6b100d4e76f49e3d2.html#other_visualizations 

#Quick lib ----

library("devtools")
library(phyloseq)
library(microbiome)

#Load data ----

pseq <- Marissa_mb2021_filtered_20240203
pseq <-  Marissa_MU42022_rare

pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))

#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#create objects

pseq.rel <- microbiome::transform(pseq, "compositional")

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
pseq.rel <- subset_samples(pseq.rel, !Age %in% c("3 dpf", "1 dpf", "18 dpf"))

#Sort the OTUs by abundance and pick the top 20
top10OTU.names = names(sort(taxa_sums(pseq.rel), TRUE)[1:10])

#Cut down the physeq.tree data to only the top 10 Phyla
top10OTU = prune_taxa(top10OTU.names, pseq.rel)

top30OTU

plot_heatmap(top60OTU)

plot_heatmap(top60OTU, sample.label="Treatment", sample.order="Age", taxa.label="Family", taxa.order="Order", low="grey", high="red", na.value="black", facet_wrap(~Age))

plot_heatmap(top10OTU, sample.label = "Treatment", sample.order = "Treatment", taxa.label = "Family", taxa.order = "Order", 
             low = "deepskyblue2", high = "red", na.value = "black") + facet_grid(~Age)


plot_heatmap(top20OTU, sample.label="Age", sample.order="Age", taxa.label="Phylum", taxa.order="Family", low="white", high="purple", na.value="grey")

plot_heatmap(top20OTU, "NMDS", "bray", title="Bray-Curtis")


#Heat map of specific bacterial groups - starting with relative data

gpac <- subset_taxa(pseq.rel, Family=="Vibrionaceae")
gpac2 <- subset_samples(gpac, Age=="Spat")
View(gpac)
gpac2 <- psmelt(gpac)

#remove outlier samples
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

ggplot(gpac2, aes(Sample, OTU, fill= Abundance)) + 
  geom_tile() +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1))

ggplot(gpac, aes(x = Sample, y = OTU, fill = Abundance)) + 
  geom_tile() +
  scale_x_discrete(limits = unique(gpac$Sample[order(gpac$Age)]), labels = gpac$Treatment) +
  xlab("Treatment") +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1))

ggplot(gpac2, aes(x = Sample, y = OTU, fill = Abundance)) + 
  geom_tile() +
  scale_x_discrete(limits = unique(gpac2$Sample[order(gpac2$Treatment)]), breaks = gpac2$Sample, labels = gpac2$Treatment) +
  xlab("Treatment") +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
  ggtitle("Vibrionaceae - Spat")

ggplot(gpac2, aes(x = Sample, y = OTU, fill = Abundance)) + 
  geom_tile() +
  scale_x_discrete(limits = unique(gpac2$Sample[order(gpac2$Treatment)]), breaks = gpac2$Sample, labels = gpac2$Sample) +
  xlab("Treatment") +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
  ggtitle("Vibrionaceae - Spat")

ggplot(gpac2, aes(x = Sample, y = OTU, fill = Abundance)) + 
  geom_tile() +
  xlab("Treatment") +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
  ggtitle("Vibrionaceae - Spat")


ggplot(gpac2, aes(x = Sample, y = OTU, fill = Abundance)) + 
  geom_tile() +
  facet_grid(~Age)
  scale_x_discrete(limits = unique(gpac2$Sample[order(gpac2$Age)]), breaks = gpac2$Sample, labels = gpac2$Treatment) +
  xlab("Treatment") +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
  ggtitle("Vibrionaceae")
  
  ggplot(gpac2, aes(x = Sample, y = OTU, fill = Abundance)) + 
    geom_tile() +
    xlab("Treatment") +
    scale_x_discrete(limits = unique(gpac2$Sample[order(gpac2$Age)]), breaks = gpac2$Sample, labels = gpac2$Treatment) +
    theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
    ggtitle("Vibrionaceae")
  
  View()
