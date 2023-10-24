#phyloseq network analysis
#source: https://joey711.github.io/phyloseq/plot_network-examples.html 
#2023-10-05
#By Marissa WL

#Quick lib ----

library("devtools")
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)


#Load data ----

Marissa_MU42022_rare <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rare.rds")

pseq <-  Marissa_MU42022_rare

#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#remove NA - here it is algae samples

pseq = subset_samples(pseq, !is.na(Treatment))

#check na removed

View(data@sam_data)

#want to remove other samples? e.g., look at:

#day 1 only

pseq <- subset_samples(pseq, !Age %in% c("Spat", "Day 03", "Day 06", "Day 15"))

#day 3 only

pseq <- subset_samples(pseq, !Age %in% c("Spat", "Day 01", "Day 06", "Day 15"))

#day 6 only

pseq <- subset_samples(pseq, !Age %in% c("Spat", "Day 03", "Day 01", "Day 15"))

#day 15 only

pseq <- subset_samples(pseq, !Age %in% c("Spat", "Day 03", "Day 06", "Day 01"))

#spat only

pseq <- subset_samples(pseq, !Age %in% c("Day 01", "Day 03", "Day 06", "Day 15"))


#plotting ----

data <- pseq

set.seed(711L)

plot_net(data, maxdist = 0.4, point_label = "Sample.ID")

plot_net(data, maxdist = 0.7, color = "Treatment", shape="Age")

#plot network function


ig <- make_network(data, max.dist=0.7)

plot_network(ig, data)

plots <- plot_network(ig, data, color="Treatment", line_weight=0.4, label=NULL)


plots


#replace automatic Jaccard distance with Brays-curtis

ig <- make_network(data, dist.fun="bray", max.dist=0.7)

plots <- plot_network(ig, data, color="Treatment", line_weight=0.4, label=NULL)


#Arrange plots as a grid ----

library(ggplot2)
library(gridExtra)

# Create your individual ggplot plots
plot1 <- ggplot(ggtitle("Day 1"))
plot3 <- ggplot(ggtitle("Day 3"))
plot6 <- ggplot(ggtitle("Day 6"))
plot15 <- ggplot(ggtitle("Day 15"))
plotas <- ggplot(ggtitle("Spat"))


# Arrange the plots in a grid
grid.arrange(plot1, plot3, plot6, plot15, plotas, ncol = 2)  
# ncol specifies the number of columns in the grid