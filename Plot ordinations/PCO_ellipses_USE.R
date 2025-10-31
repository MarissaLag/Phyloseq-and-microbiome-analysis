#PCO plots faceted by age adn with ellipses
#note: ellipses are not confidence intervals (have no statistical value)

## Load libraries
library(devtools)
library(phyloseq)
library(microbiome)
library(ggalt)
library(ggplot2)
library(gridExtra)
library(ggpubr)

#set theme
theme.marissa <- function() {
  theme_classic(base_size = 14) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"))
}

theme_set(theme.marissa())


#load RDS and filter if needed ----

#MU42022 filtering
pseq <- MU42022_filtered_Oct92024
#pseq <- MU42022_filtered_NOT_rarefied
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")

#MB2021 filtering
pseq <- Marissa_mb2021_filtered_20240203
pseq <- mb2021_filtered_NOT_rarefied_normalized
pseq <- subset_samples(pseq, Age %in% c("18 dpf"))
pseq <- subset_samples(pseq, !Family %in% c("9"))

#Mu42022 - combining to add algae
merged_phyloseq <- merge_phyloseq(MU42022_filtered_algae, MU42022_filtered_Oct92024)
View(merged_phyloseq@otu_table)
pseq <- merged_phyloseq

#PB2023
pseq <- PB2023_spat_10X_limited_CSS
pseq <- PB2023_rarefied_3000
pseq <- PB2023_rarefied_1313
pseq <- subset_samples(pseq, !Treatment %in% c("Continuous Probiotics", "James"))

#Sam
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Age %in% c("18 dpf"))

#convert to compositional data

pseq.rel <- microbiome::transform(pseq, "compositional")

#plot MDS/PcoA ----
#Create PCO for each time-point then combine plots **must be done for correct ordination
#note: different projects have different time-points
#note: if a treatment is lost at a certain time-point, adjust colours manually 

set.seed(4235421)

#info on geom_encircle
?ggalt::geom_encircle

#MU42022 ----
#All time-points
pseq <- Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% "Spat")
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

p_legend <- plot_ordination(pseq.rel, ord, color = "Microbial.Source", shape = "Day") +
  #geom_point(aes(fill = Treatment), size = 6) +
  geom_point(size = 8) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen", "yellow")) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen", "yellow")) + # Matching colors for ellipses and points
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  #geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)

p_legend

p <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) + # Matching colors for ellipses and points
  ggtitle("All time-points") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  #geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
p

#To add circles around time-points
p <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) +
  ggtitle("All time-points") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_encircle(aes(group = Age), expand = 0.2, color = "black", alpha = 0.5, fill = NA)  # Removed fill aesthetic mapping

# Show the plot
print(p)

#1 dpf
#important to note that the PB sequence was not removed
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
#pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature")) #Add/remove "orange" if HT included
pseq <- subset_samples(pseq, Age %in% c("Day 01"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p1 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Family", label = "Sample.ID") +
  geom_point(aes(fill = Treatment), size = 9) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + # Matching colors for ellipses and points
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
  #geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
p1


p1 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Family") +
  geom_point(aes(fill = Treatment), size = 8) +
  geom_text(aes(label = Sample.ID.seq), size = 5, vjust = -1, check_overlap = TRUE) +  # Add this line to label points
  scale_colour_manual(values = c("darkgrey", "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

print(p1)


#3 dpf
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
#pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))
pseq <- subset_samples(pseq, Age %in% c("Day 03"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p3 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), size = 6, shape = 24) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + # Matching colors for ellipses and points
  ggtitle("Day 3") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
#geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
p3


#6 dpf
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
#pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))
pseq <- subset_samples(pseq, Age %in% c("Day 06"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p6 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), size = 6, shape = 22) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey",  "cornflowerblue", "orange")) + # Matching colors for ellipses and points
  ggtitle("Day 6") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p6

#15 dpf
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))
pseq <- subset_samples(pseq, Age %in% c("Day 15"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p15 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), size = 6, shape = 3) +
  scale_colour_manual(values = c("darkgrey", "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + # Matching colors for ellipses and points
  ggtitle("Day 15") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p15

#Spat
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
#pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% c("Spat"))
pseq@sam_data$Genetics <- as.factor(pseq@sam_data$Genetics)
pseq.rel <- microbiome::transform(pseq, "compositional")

#note: for MU42022 no Heat treatment by Spat (blue colour)

ord <- ordinate(pseq.rel, "MDS", "bray")

pSpat <- plot_ordination(pseq.rel, ord, color = "Treatment") +
  geom_point(aes(fill = Treatment), size = 6, shape = 7) +
  scale_colour_manual(values = c("darkgrey", "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + # Matching colors for ellipses and points
  ggtitle("Spat") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

pSpat

pSpat <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Genetics", label = "Genetics") +
  geom_point(aes(fill = Treatment), size = 6) +
  scale_colour_manual(values = c("darkgrey", "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + # Matching colors for ellipses and points
  ggtitle("Spat") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")



#arrange PCO plots -----

grid.arrange(p1, p3, p6, p15, pSpat, ncol = 3)

ggarrange(p1, p3, p6, p15, pSpat, nrow = 2, common.legend = TRUE, legend="bottom")

# Arrange all plots, removing legends from all except the last plot (pSpat)
ggarrange(p1 + theme(legend.position = "none"), 
          p3 + theme(legend.position = "none"),
          p6 + theme(legend.position = "none"),
          p15 + theme(legend.position = "none"),
          pSpat,  # Use pSpat as the source of the common legend
          ncol = 3,  # Set 3 columns to make sure all plots fit
          nrow = 2,  # Specify 2 rows to make the layout balanced
          common.legend = TRUE, 
          legend = "bottom")


#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  
  if (length(leg) == 0) {
    stop("Legend not found in the plot.")
  }
  
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend<-g_legend(p_legend)

#adjust legend
mylegend <- mylegend + theme(legend.text = element_text(size = 8))



p_combined <- grid.arrange(
  arrangeGrob(
    p + theme(legend.position = "none"),
    p1 + theme(legend.position = "none"),
    p3 + theme(legend.position = "none"),
    p6 + theme(legend.position = "none"),
    p15 + theme(legend.position = "none"),
    pSpat + theme(legend.position = "none"),
    ncol = 3
  ),
  mylegend,
  ncol = 1,
  heights = c(4, 1)  # Adjust the height ratio between plots and legend
)

# Show the combined plot
print(p_combined)

#PB2023 ----
pseq <- PB2023_spat_10X_limited_CSS
pseq <- subset_samples(pseq, !Treatment %in% c("Continuous Probiotics", "James"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

pSpat <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Family") +
  geom_point(aes(fill = Treatment), size = 8) +
  scale_colour_manual(values = c("darkgrey", "orange", "forestgreen")) +
  scale_fill_manual(values = c("darkgrey", "orange", "forestgreen")) + # Matching colors for ellipses and points
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

pSpat

#mb2021 ----

#Trying different ordination methods and distance methods
#Above code uses MDS/PCOa plots with Brays-Curtis distance matrix (relative abundance)


#NMDS = Performs Non-metric MultiDimenstional Scaling of a sample-wise ecological distance matrix 
#onto a user-specified number of axes, k.

pseq <- mb2021_filtered_NOT_rarefied_normalized
#pseq <- Marissa_mb2021_filtered_20240203
pseq <- subset_samples(pseq, Age %in% c("Spat"))
pseq <- subset_samples(pseq, !Family %in% c("9"))

#pseq <- subset_samples(pseq, !Family %in% c("1"))

#correct family column for mb2021

pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")


p <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Family") +
  geom_point(aes(fill = Treatment), size = 6) +
  ggforce::geom_mark_ellipse(aes(color = Treatment)) +
  scale_colour_manual(values = c("#F8766D", "chartreuse3", "deepskyblue3")) +
  ggtitle("All time-points") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",
        axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold")) +
  xlim(-0.6, 0.4) +
  ylim(-0.4, 0.35)

p

p1 <- plot_ordination(pseq.rel, ord, color = "Treatment") +
  geom_point(aes(fill = Treatment), shape = 16, size = 6) +
  ggforce::geom_mark_ellipse(aes(color = Treatment)) +
  scale_colour_manual(values = c("#F8766D", "chartreuse3", "deepskyblue3")) +
  #ggtitle("1 dpf") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
        axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold")) +
  xlim(-0.35, 0.4) +
  ylim(-0.4, 0.3) +
  facet_wrap(~Age)
p1


p2 <- plot_ordination(pseq.rel, ord, color = "Treatment") +
  geom_point(aes(fill = Treatment), shape = 17, size = 6) +
  ggforce::geom_mark_ellipse(aes(color = Treatment)) +
  scale_colour_manual(values = c("#F8766D", "chartreuse3", "deepskyblue3")) +
  #ggtitle("18 dpf") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
        axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold")) +
  ylim(-0.3, 0.55) +
  xlim(-0.5, 0.75) +
  facet_wrap(~Age)
p2

p3 <- plot_ordination(pseq.rel, ord, color = "Treatment") +
  geom_point(aes(fill = Treatment), shape = 15, size = 6) +
  ggforce::geom_mark_ellipse(aes(color = Treatment)) +
  scale_colour_manual(values = c("#F8766D", "chartreuse3", "deepskyblue3")) + 
  #ggtitle("Spat") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
        axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold")) +
  ylim(-0.3, 0.3) +
  xlim(-0.45, 0.3) +
  facet_wrap(~Age) 
  
p3


ggarrange(p1, p2, p3, nrow = 1, ncol =3)

ggarrange(p1, p2, p3, nrow = 1, ncol =3, common.legend = TRUE, legend="bottom")

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  
  if (length(leg) == 0) {
    stop("Legend not found in the plot.")
  }
  
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend<-g_legend(p_legend)

#adjust legend
mylegend <- mylegend + theme(legend.text = element_text(size = 8))



p_combined <- grid.arrange(
  arrangeGrob(
    p + theme(legend.position = "none"),
    p1 + theme(legend.position = "none"),
    p3 + theme(legend.position = "none"),
    p6 + theme(legend.position = "none"),
    p15 + theme(legend.position = "none"),
    pSpat + theme(legend.position = "none"),
    ncol = 3
  ),
  mylegend,
  ncol = 1,
  heights = c(4, 1)  # Adjust the height ratio between plots and legend
)

# Show the combined plot
print(p_combined)


#Sam ----

#3 dpf
#important to note that the PB sequence was not removed
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Day %in% c("3"))
pseq <- subset_samples(pseq, Organism %in% c("Oyster"))
pseq <- subset_samples(pseq, !Sample.ID %in% c("20r3", "20r2"))
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

library(ggrepel)

p1 <- plot_ordination(pseq.rel, ord, color = "Microbial.Source", shape = "Genetic.Background") +
  geom_point(aes(fill = Microbial.Source), size = 8) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + 
  #geom_text_repel(aes(label = Sample.ID), size = 3) +  # add labels with repelling
  ggtitle("Water samples") +
  theme(plot.title = element_text(hjust = 0.5))
p1

#6 dpf #only water samples
#important to note that the PB sequence was not removed
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Day %in% c("6"))
#pseq <- subset_samples(pseq, Organism %in% c("Oyster"))
#pseq <- subset_samples(pseq, !Sample.ID %in% c("20r3", "20r2"))
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

p2 <- plot_ordination(pseq.rel, ord, color = "Microbial.Source", shape = "Organism") +
  geom_point(aes(fill = Microbial.Source), size = 8) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + 
  ggtitle("Day 6") +
  theme(plot.title = element_text(hjust = 0.5))

p2

#10 dpf #only larval samples
#important to note that the PB sequence was not removed
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Day %in% c("10"))
#pseq <- subset_samples(pseq, Organism %in% c("Oyster"))
#pseq <- subset_samples(pseq, !Sample.ID %in% c("20r3", "20r2"))
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

p3 <- plot_ordination(pseq.rel, ord, color = "Microbial.Source", shape = "Genetic.Background") +
  geom_point(aes(fill = Microbial.Source), size = 8) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + 
  ggtitle("Day 10") +
  theme(plot.title = element_text(hjust = 0.5))

p3

#15 dpf #only larval samples
#important to note that the PB sequence was not removed
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Day %in% c("15"))
#pseq <- subset_samples(pseq, Organism %in% c("Oyster"))
#pseq <- subset_samples(pseq, !Sample.ID %in% c("20r3", "20r2"))
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

p4 <- plot_ordination(pseq.rel, ord, color = "Microbial.Source", shape = "Genetic.Background") +
  geom_point(aes(fill = Microbial.Source), size = 8) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + 
  ggtitle("Day 15") +
  theme(plot.title = element_text(hjust = 0.5))
p4

#Spat
#important to note that the PB sequence was not removed
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Day %in% c("Spat"))
#pseq <- subset_samples(pseq, Organism %in% c("Oyster"))
#pseq <- subset_samples(pseq, !Sample.ID %in% c("20r3", "20r2"))
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

p5 <- plot_ordination(pseq.rel, ord, color = "Microbial.Source", shape = "Genetic.Background") +
  geom_point(aes(fill = Microbial.Source), size = 8) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + 
  ggtitle("Spat") +
  theme(plot.title = element_text(hjust = 0.5))
p5

#Spat 2
#important to note that the PB sequence was not removed
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Day %in% c("Spat2"))
#pseq <- subset_samples(pseq, Organism %in% c("Oyster"))
#pseq <- subset_samples(pseq, !Sample.ID %in% c("20r3", "20r2"))
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

p6 <- plot_ordination(pseq.rel, ord, color = "Microbial.Source", shape = "Genetic.Background") +
  geom_point(aes(fill = Microbial.Source), size = 8) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + 
  ggtitle("Spat 2") +
  theme(plot.title = element_text(hjust = 0.5))
p6


#Water samples
pseq <- Sam_all_samples_partial_rare_CSS
pseq <- subset_samples(pseq, Organism %in% c("Water"))
#pseq <- subset_samples(pseq, !Sample.ID %in% c("20r3", "20r2"))
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

p6 <- plot_ordination(pseq.rel, ord, color = "Microbial.Source", shape = "Day") +
  geom_point(aes(fill = Microbial.Source), size = 8) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) + 
  ggtitle("Water") +
  theme(plot.title = element_text(hjust = 0.5))
p6


#Combine plots
ggarrange(p1, p3,p4, p5, p6, nrow = 3, ncol =2, common.legend = TRUE, legend="bottom")


