#Permanova & PCO plot

##Quick lib ---- 
#(if you have packages downloaded already)

library("devtools")
library(phyloseq)
library(microbiome)
library(vegan)

#Load and name data ----

Marissa_MU42022_rarefied_20231016 <- readRDS("~/GitHub/Phyloseq and microbiome analysis/Marissa_MU42022_rarefied_20231016.rds")

pseq <- Marissa_MU42022_rarefied_20231016

pseq <- Marissa_mb2021_filtered_20240203

#for mb2021 remove 3 dpf

pseq <- subset_samples(pseq, !Age %in% "3 dpf")

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

summary <- adonis2(pseq_bray ~ Treatment*Age*Family, data = metadata)


#Homogeneity of dispersion test
beta <- betadisper(pseq_bray, metadata$Family)
permutest(beta)

##treatment/family follows homogeneity but Age does not (p = 0.001)


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







