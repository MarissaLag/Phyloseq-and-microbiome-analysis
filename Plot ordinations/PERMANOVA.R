#Permanova & PCO plot

##Quick lib ---- 
#(if you have packages downloaded already)

library("devtools")
library(phyloseq)
library(microbiome)
library(vegan)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library(pairwiseAdonis)

#Load and name data ----

Marissa_MU42022_rarefied_20231016 <- readRDS("~/GitHub/Phyloseq and microbiome analysis/Marissa_MU42022_rarefied_20231016.rds")

#pseq <- Marissa_MU42022_rarefied_20231016

pseq <- Marissa_mb2021_filtered_20240203

pseq <- mb2021_filtered_NOT_rarefied

#for mb2021 remove 3 dpf and T9

pseq <- subset_samples(pseq, !Family %in% "9")

pseq <- subset_samples(pseq, Age %in% "Spat")

pseq <- subset_samples(pseq, !Library_Name %in% c("T9r1", "T9r3"))

#compositional

pseq <- microbiome::transform(pseq, "compositional")

#correct family column

# Replace values according to your mapping
pseq@sam_data$Family <- ifelse(pseq@sam_data$Family %in% c(9, 13), "1",
                               ifelse(pseq@sam_data$Family %in% c(10, 14), "2",
                                      ifelse(pseq@sam_data$Family %in% c(11, 15), "3",
                                             ifelse(pseq@sam_data$Family %in% c(12, 16), "4", pseq@sam_data$Family))))

# Print the updated Family column
print(pseq@sam_data$Family)
View(pseq@sam_data) 


#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#PERMAOVA
#source: https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
set.seed(452)

pseq_bray <- phyloseq::distance(pseq, method = "bray")

metadata <- as(sample_data(pseq), "data.frame")

summary <- adonis2(pseq_bray ~ Treatment*Family, data = metadata)

summary


#Homogeneity of dispersion test
beta <- betadisper(pseq_bray, metadata$Family)
permutest(beta)

##treatment/family follows homogeneity but Age does not (p = 0.001)

#Post-Hoc

pairwise.adonis(pseq_bray, phyloseq::sample_data(pseq)$Treatment)

#CAP plots ----

pseq <- Marissa_mb2021_filtered_20240203

#for mb2021 remove 3 dpf

pseq <- subset_samples(pseq, !Age %in% "3 dpf")
pseq_bray <- phyloseq::distance(pseq, method = "bray")

# CAP ordinate
cap_ord <- ordinate(
  physeq = pseq, 
  method = "CAP",
  distance = pseq_bray,
  formula = ~Age 
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = pseq, 
  ordination = cap_ord, 
  color = "Treatment", 
  axes = c(1,2)
) + 
  aes(shape = Age) + 
  geom_point(aes(colour = Treatment), alpha = 0.4, size = 4) + 
  #geom_point(colour = "grey90", size = 3) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
                                "#1919ff", "darkorchid3", "magenta")
  )


arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrowdf$labels <- ifelse(grepl("High", arrowdf$labels), "High salinity", "Low salinity")
print(arrowdf)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.01, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "grey", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE,
    nudge_x = 0.05,
    nudge_y = 0.05
  ) +
  facet_wrap(~Age)

##need to find code for pairwise comparison




#PERMANOVA of specific microbiota #Code below doesn't work...can't make brays similarity matrix with only Vibrio?

pseq <- Marissa_mb2021_filtered_20240203

#for mb2021 remove 3 dpf and T9

pseq <- subset_samples(pseq, !Age %in% "3 dpf")

pseq <- subset_samples(pseq, Age %in% "18 dpf")

pseq <- subset_samples(pseq, !Library_Name %in% c("T9r1", "T9r3"))

pseq <- microbiome::transform(pseq, "compositional")


subset <- subset_taxa(pseq, Family=="Vibrionaceae")
View(subset@otu_table)

set.seed(452)

pseq_bray <- phyloseq::distance(subset, method = "bray")

metadata <- as(sample_data(subset), "data.frame")

summary <- adonis2(pseq_bray ~ Treatment*Age, data = metadata)

summary


#SIMPER

library(phyloseq)
library(vegan)

pseq<- Marissa_mb2021_filtered_20240203
#filter data (if needed)
pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))
pseq <- subset_samples(pseq, Age %in% c("1 dpf"))

pseq <- psmelt(pseq)

# Define groupings for SIMPER analysis (e.g., Treatment)
grouping <- pseq@sam_data$Treatment

otu <- pseq@otu_table

# Perform SIMPER analysis
simper_result <- simper(otu, grouping)

# Print the results
print(simper_result)

#test each ASV for significance

# Extract OTU abundance matrix
otu_matrix <- as.matrix(otu_table(pseq))

# Transpose the OTU matrix so that ASVs are in rows and samples are in columns
otu_matrix <- t(otu_matrix)

# Perform Kruskal-Wallis test for each ASV
kruskal_results <- apply(otu_matrix, 1, function(x) kruskal.test(x ~ grouping))

# Extract p-values from the test results
p_values <- sapply(kruskal_results, function(x) x$p.value)

# Adjust p-values for multiple testing if needed (e.g., using Bonferroni correction)

adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Print the results
results <- data.frame(ASV = rownames(otu_matrix), p_value = adjusted_p_values)
print(results)


#different post hoc test
#source:https://www.yanh.org/2021/01/01/microbiome-r/ 

metadata <- data.frame(sample_data(ps.rarefied))
test.adonis <- adonis(dist ~ body.site, data = metadata)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis

cbn <- combn(x=unique(metadata$body.site), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps.rarefied, body.site %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis(phyloseq::distance(ps.subs, method = "bray") ~ body.site, 
                               data = metadata_sub)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

