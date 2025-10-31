#Permanova & PCO plot

##Quick lib ---- 
#(if you have packages downloaded already)

library("devtools")
library(phyloseq)
# library(BiocManager)
# BiocManager::install("microbiome")
library(microbiome)
library(vegan)
library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", force = TRUE) 
library(pairwiseAdonis)

#Load and name data ----

pseq <- MU42022_filtered_Oct92024

pseq <- Sam_all_samples_partial_rare_CSS

#MU42022 - remove F4 and algal sample

pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

# pseq <- Marissa_mb2021_filtered_20240203
# 
# pseq <- mb2021_filtered_NOT_rarefied
# 
# pseq <- mb2021_filtered_NOT_rarefied_normalized

#for mb2021 remove 3 dpf and T9

# pseq <- subset_samples(pseq, !Family %in% "9")
# 
# pseq <- subset_samples(pseq, Age %in% "Day 01")
# 
# pseq <- subset_samples(pseq, !Library_Name %in% c("T9r1", "T9r3"))

#compositional

pseq <- microbiome::transform(pseq, "compositional")

#look at data

p <- plot_landscape(pseq, method = "NMDS", distance = "bray", col = "Treatment", size = 3)
print(p)


#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#make treatment a factor

Metadata$Genetic.Background <- as.factor(Metadata$Genetic.Background)
Metadata$Microbial.Source <- as.factor(Metadata$Microbial.Source)
Metadata$Organism <- as.factor(Metadata$Organism)
Metadata$Day <- as.factor(Metadata$Day)

#PERMAOVA
#source: https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Stage %in% c("Parent"))
pseq <- subset_samples(pseq, Organism %in% c("Oyster"))
set.seed(452)

pseq_bray <- phyloseq::distance(pseq, method = "bray", weighted = TRUE)

Metadata <- as(sample_data(pseq), "data.frame")

Metadata$Genetic.Background <- as.factor(Metadata$Genetic.Background)
Metadata$Microbial.Source <- as.factor(Metadata$Microbial.Source)

summary_additive <- adonis2(pseq_bray ~ Microbial.Source + Genetic.Background,
                            data = Metadata, by = "margin")
summary_additive

#Homogeneity of dispersion test
beta <- betadisper(pseq_bray, Metadata$Microbial.Source)
permutest(beta)

##treatment/family follows homogeneity but Age does not (p = 0.001)

#Post-Hoc
#Issues with pairwise.adonis - not enough permutations

pairwise.adonis(pseq_bray, Metadata$Microbial.Source)

pairwise.adonis(pseq_bray, Metadata$Treatment)

#Cannot do Tukey's post hoc test on permanova (must be anova)

#Because odd pairwise tests (apparently because not there are not enough permutations to run a comparison)
#Trying a Negative Binomial generalized linear model or MANOVA (assumes data is normal)



#CAP plots ----
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Day %in% c("Spat2"))
pseq <- subset_samples(pseq, Organism %in% c("Oyster"))
pseq_bray <- phyloseq::distance(pseq, method = "bray", weighted = TRUE)

# CAP ordinate
cap_ord <- ordinate(
  physeq = pseq, 
  method = "CAP",
  distance = pseq_bray,
  formula = ~Microbial.Source
)

# CAP plot
p5 <- plot_ordination(
  physeq = pseq,
  ordination = cap_ord,
  color = "Microbial.Source",
  axes = c(1,2)
) +
  aes(shape = Genetic.Background) +
  geom_point(aes(colour = Microbial.Source),size = 8) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) +
  ggtitle("Spat2 - CAP plot") +
  theme(plot.title = element_text(hjust = 0.5))
p5

#Combine plots
ggarrange(p1, p2,p3, p4, p5, nrow = 2, ncol =3, common.legend = TRUE, legend="bottom")



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

