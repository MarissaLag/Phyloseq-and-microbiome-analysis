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

pseq <- mb2021_filtered_NOT_rarefied_normalized

#for mb2021 remove 3 dpf and T9

pseq <- subset_samples(pseq, !Family %in% "9")

pseq <- subset_samples(pseq, Age %in% "1 dpf")

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



#look at data

p <- plot_landscape(pseq, method = "NMDS", distance = "bray", col = "Treatment", size = 3)
print(p)

#remove ASVs that have zero abundance in your subset of data

otu_table <- phyloseq::otu_table(pseq)
otu_matrix <- as(otu_table, "matrix")
taxa_to_keep <- rowSums(otu_matrix) > 0
otu_table_filtered <- otu_table[taxa_to_keep, ]
pseq_filtered <- phyloseq::otu_table(otu_table_filtered, taxa_are_rows = TRUE)
pseq <- merge_phyloseq(pseq, pseq_filtered)


#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#make treatment a factor

Metadata$Treatment <- as.character(Metadata$Treatment)

#PERMAOVA
#source: https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
set.seed(452)

pseq_bray <- phyloseq::distance(pseq, method = "bray", weighted = TRUE)

metadata <- as(sample_data(pseq), "data.frame")

summary <- adonis2(pseq_bray ~ Treatment*Family, data = metadata)

summary




#Homogeneity of dispersion test
beta <- betadisper(pseq_bray, metadata$Treatment)
permutest(beta)

##treatment/family follows homogeneity but Age does not (p = 0.001)

#Post-Hoc

pairwise.adonis(pseq_bray, phyloseq::sample_data(pseq)$Family)

pairwise.adonis(pseq_bray, Metadata$Treatment)

#getting weird results for some of the pairwise tests....trying Tukeys instead

pseq_bray <- phyloseq::distance(pseq, method = "bray")

# Extract metadata and create a data frame
metadata <- as(sample_data(pseq), "data.frame")

# Calculate beta diversity dispersion
bd <- betadisper(pseq_bray, metadata$Treatment)

# Perform ANOVA using aov
bd_aov <- aov(bd$distances ~ metadata$Treatment) 

# Perform Tukey's HSD post-hoc test
tukey_result <- TukeyHSD(bd_aov)

# View results
print(tukey_result)


#Because odd pairwise tests (apparently because not there are not enough permutations to run a comparison)
#Trying a Negative Binomial generalized linear model or MANOVA (assumes data is normal)






#CAP plots ----

pseq <- Marissa_mb2021_filtered_20240203


# CAP ordinate
cap_ord <- ordinate(
  physeq = pseq, 
  method = "CAP",
  distance = pseq_bray,
  formula = ~Treatment
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = pseq, 
  ordination = cap_ord, 
  color = "Treatment", 
  axes = c(1,2)
) + 
  aes(shape = Treatment) + 
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

