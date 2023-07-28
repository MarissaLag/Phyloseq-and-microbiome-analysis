##mb2021 project: transfer Qiime data into R -> coordinates -> PCoA and MDS


install.packages("devtools")

library("devtools")

#install source: https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDbData")

BiocManager::install("phyloseq")
##update all changes
#if phyloseq not installing use below code
BiocManager::install("phyloseq", force = TRUE)

library(phyloseq)

#source: https://www.bioconductor.org/packages/release/bioc/html/microbiome.html
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")
##update all changes
library(microbiome)

#load data

Marissa_Osyter <- readRDS("~/mb2021/Marissa_Osyter.rds")

Marissa_Osyter

#create objects

OTU = Marissa_Osyter@otu_table
Tax = Marissa_Osyter@tax_table
Metadata = Marissa_Osyter@sam_data
Tree = Marissa_Osyter@phy_tree

##rarify seqs to exclude seqs that contribute to less than 0.05% of abundance data

physeq <- phyloseq(OTU)

sequence_abundance <- taxa_sums(physeq)

threshold <- 0.0005

# Get the names of sequences that meet the threshold
abundant_sequences <- names(sequence_abundance[sequence_abundance > threshold])

# Filter out sequences not meeting the threshold and create a new phyloseq object
physeq_rarified <- prune_taxa(abundant_sequences, physeq)

##rarify data to make sequences an even depth

Rare_OTU <-rarefy_even_depth(physeq_rarified, sample.size= 5000)

##12 samples removedbecause they contained fewer reads than `sample.size`.
##Up to first five removed samples are: 
  
#F1C3F1H3F2H3F2L1F2L3 (i.e., most of day 3)
...
#6397OTUs were removed because they are no longer 
#present in any sample after random subsampling


##Only one sample from day 3 left - remove it from analysis 

excluded_samples <- "F4L3"

Rare_OTU_filtered <- Rare_OTU[, !colnames(Rare_OTU) %in% excluded_samples]

# Convert the rarefied OTU matrix to a new phyloseq object and replace old OTU table in Marissa_Oyster data - Yay! I can tell it worked because sample number went from 63 to 50 (12 + 1 samples removed above)

physeq_rarified_OTU <- phyloseq(otu_table(Rare_OTU_filtered))

Marissa_Oyster@otu_table <- physeq_rarified_OTU


# Save the modified marissa_oyster_data back to the RDS file

saveRDS(Marissa_Oyster, "C:/Users/maris/OneDrive/Documents/USRA2021/mb2021/Data/Marissa_Oyster.rds")


#create objects

pseq <- Marissa_Osyter

pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq, level = "Genus", detection = 50/100, prevalence = 70/100)

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

#PERMANOVA - source: https://microbiome.github.io/tutorials/PERMANOVA.html

library(vegan)
meta <- meta(pseq.rel)
permanova <- adonis2(t(OTU) ~ Tank_treatment,
                     data = meta, permutations=999, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["Tank_treatment", "Pr(>F)"])
##getting result "NULL"....? not working

##carrying on with compositional analysis

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq, ord, color = "Tank_treatment", shape = "Age") + geom_point(size = 4)

##change colour palette for PCOA plot

names(Data)[names(Data) == "treatment"] <- "Treatment"

custom_color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00") # Example color palette

plot_ordination(pseq, ord, color = "Tank_treatment", shape = "Age") +
  geom_point(size = 4) +
  scale_color_manual(values = custom_color_palette, labels = c("Control", "High salinity", "Low salinity")) +
  scale_shape_manual(values = c("day_1" = 16, "day_18" = 17, "day_3" = 15, "spat" = 3), 
                     labels = c("Day 1", "Day 18", "Day 3", "Spat"))

##exlude day 3 data - issue = "NA" appears under tank_treatment legend

filtered_data <- subset_samples(pseq, !grepl("day_3", Age))

# Use plot_ordination with the filtered data
plot_ordination(filtered_data, ord, color = "Tank_treatment", shape = "Age") +
  geom_point(size = 4) +
  scale_color_manual(values = custom_color_palette, labels = c("Control", "High salinity", "Low salinity")) +
  scale_shape_manual(values = c("day_1" = 16, "day_18" = 17, "spat" = 3), 
                     labels = c("Day 1", "Day 18", "Spat"))


##Conical correspondance analysis
##source; https://microbiome.github.io/tutorials/Ordination.html
Kindly cite this work as follows: "Leo Lahti, Sudarshan Shetty et al. (2017). Tools for microbiome analysis in R. Version . URL: http://microbiome.github.com/microbiome. Check also the relevant references listed in the manual page of each function.

The package utilizes tools from a number of other R extensions, including dplyr (Wickham, Francois, Henry, and MÃ¼ller, 2017), ggplot2 (Wickham, 2009), phyloseq (McMurdie and Holmes, 2013), tidyr (Wickham, 2017), vegan (Oksanen, Blanchet, Friendly, Kindt, Legendre, McGlinn, Minchin, O’Hara, Simpson, Solymos, Stevens, Szoecs, and Wagner, 2017).

# With samples
pseq.cca <- ordinate(pseq, "CCA")
p <- plot_ordination(pseq, pseq.cca,
                     type = "samples", color = "Tank_treatment", shape = "Age")
p <- p + geom_point(size = 4)
print(p)

# With taxa:
p <- plot_ordination(pseq, pseq.cca,
                     type = "taxa", color = "Phylum")
p <- p + geom_point(size = 4)
print(p)

#split plot

plot_ordination(pseq, pseq.cca,
                type = "split", shape = "Tank_treatment", 
                color = "Phylum", label = "Tank_treatment")

#t-SNE

library(vegan)
library(microbiome)
library(Rtsne)

set.seed(423542)

method <- "tsne"
trans <- "hellinger"
distance <- "euclidean"

# Distance matrix for samples
ps <- microbiome::transform(pseq, trans)

# Calculate sample similarities
dm <- vegdist(otu_table(ps), distance)

# Run TSNE
tsne_out <- Rtsne(dm, dims = 2) 
proj <- tsne_out$Y
rownames(proj) <- rownames(otu_table(ps))

library(ggplot2)
p <- plot_landscape(proj, legend = T, size = 1) 
print(p)

#alpha diversity #within sample-type variations
library(microbiome)
library(knitr)


tab <-microbiome::alpha(pseq, index = "all")
kable(head(tab))

#Richness
tab <- richness(pseq)
kable(head(tab))

#dominance (most dominant species)
tab <- dominance(pseq, index = "all")
kable(head(tab))

#rarity and low abundance
tab <- rarity(pseq, index = "all")
kable(head(tab))

#core abundance
tab <- core_abundance(pseq, detection = .1/100, prevalence = 50/100)
kable(head(tab))

#Gini index = Gini index is a common measure for inequality in economical income. The inverse gini index (1/x) can also be used as a community diversity measure.

tab <- inequality(pseq)
kable(head(tab))

#visualization - shannon index (alpha diversity) by tank treatment

p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Tank_treatment")
p.shannon <- p.shannon +
    theme_minimal() +
    labs(x = "", y = "Shannon diversity") +
    theme(
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)
    ) +
    scale_x_discrete(labels = c("Control", "High salinity", "Low salinity"))                          

p.shannon <- p.shannon +
    theme_minimal() +
    labs(x = "", y = "Shannon diversity") +
    theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold")  # Modify legend title properties
    ) +
    guides(fill = guide_legend(title = "Seawater treatment"))

p.shannon

##make text bold

p.shannon <- p.shannon +
    theme_minimal() +
    labs(x = "", y = expression(bold("Shannon diversity"))) +
    theme(
        axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold")
    ) +
    scale_x_discrete(
        breaks = c("Control", "High salinity", "Low salinity"),
        labels = c(expression(bold("Control")),
                   expression(bold("High salinity")),
                   expression(bold("Low salinity")))
    )

#visualization - shannon index (alpha diversity) by Age

p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Age")


p.shannon <- p.shannon +
    theme_minimal() +
    labs(x = "", y = "Shannon diversity") +
    theme(
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)
    ) +
    scale_x_discrete(labels = c("Day 1", "Day 18", "Day 3", "Spat"))
    
p.shannon

#testing for differences in alpha diversity by tank treatment - non-parametric Kolmogorov-Smirnov test
# Construct the data
d <- meta(pseq)
d$diversity <- microbiome::diversity(pseq, "shannon")$shannon

# Split the values by group
spl <- split(d$diversity, d$Tank_treatment)

# Kolmogorov-Smironv test
pv <- ks.test(spl$control, spl$low_salinity, spl$high_salinity)$p.value
print(pv)
# Adjust the p-value
padj <- p.adjust(pv)
print(padj)
##non-signif (p value = 0.9761322)

#testing for differences in alpha diversity by Age - non-parametric Kolmogorov-Smirnov test
# Construct the data
d <- meta(pseq)
d$diversity <- microbiome::diversity(pseq, "shannon")$shannon
# Split the values by group
spl <- split(d$diversity, d$Age)
# Kolmogorov-Smironv test
pv <- ks.test(spl$day_1, spl$day_3, spl$day_18, spl$spat)$p.value
print(pv)
# Adjust the p-value
padj <- p.adjust(pv)
print(padj)
##signif difference in alpha diversity between ages (p value = 0.002749572)

#moving on to community structure

# Make sure we use functions from correct package
transform <- microbiome::transform

# Merge rare taxa to speed up examples
pseq <- transform(pseq, "compositional")
pseq <- aggregate_rare(pseq, level = "Genus", detection = 1/100, prevalence = 50/100)

# from https://github.com/hrbrmstr/hrbrthemes
install.packages("hrbrthemes")
install.packages("gcookbook")
library(hrbrthemes)
library(gcookbook)
library(tidyverse)

p <- plot_composition(pseq.fam.rel,
                      taxonomic.level = "Phylum",
                      sample.sort = "Tank_treatment",
                      x.label = "Tank_treatment") +
  scale_fill_brewer("Phylum", palette = "Paired") +
  guides(fill = guide_legend(ncol = 4)) +
  scale_y_percent() +
  labs(x = "Tank treatment", y = "Relative abundance (%)",
       title = "Relative abundance data") + 
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
print(p) 

p <- plot_composition(pseq.fam.rel,
                      average_by = "Tank_treatment", 
                      transform = "compositional") +
  scale_fill_brewer("Genera", palette = "Paired") +
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
print(p)

#compositional heatmaps - can't figure out how to label samples on y axis

p <- microbiome::transform(pseq.fam.rel, "compositional") %>% 
  plot_composition(plot.type = "heatmap",
                   sample.sort = "Tank_treatment", 
                   otu.sort = "neatmap") +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size = 9, hjust=1),
        legend.text = element_text(size = 8)) +
  ylab("Samples")
print(p)



Rhodobacteraceae <- "D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhodobacterales;D_4__Rhodobacteraceae"


##Venn diagrams - have not figured out how to do it.

#We select for the row by name with OTU.clean[“name”,]
#We select the columns with a value >0 with OTU.clean[,apply()]

OTU.Control = colnames(OTU.clean["Control", apply(OTU.clean["Control",], MARGIN=2, function(x) any(x >0))])

OTU.High_salinity = colnames(OTU.clean["High_salinity", apply(OTU.clean["High_salinity",], MARGIN=2, function(x) any(x >0))])

OTU.Low_salinity = colnames(OTU.clean["Low_salinity",apply(OTU.clean["Low_salinity",], MARGIN=2, function(x) any(x >0))])





