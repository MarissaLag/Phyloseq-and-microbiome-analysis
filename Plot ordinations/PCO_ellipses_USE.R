#PCO plots faceted by age adn with ellipses
#note: ellipses are not confidence intervals (have no statistical value)

## Load libraries
library("devtools")
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
pseq <- Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#MB2021 filtering
pseq <- Marissa_mb2021_filtered_20240203
pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))

#convert to compositional data

pseq.rel <- microbiome::transform(pseq, "compositional")

#plot MDS/PcoA ----
#Create PCO for each time-point then combine plots **must be done for correct ordination
#note: different projects have different time-points
#note: if a treatment is lost at a certain time-point, adjust colours manually 

set.seed(4235421)

#info on geom_encircle
?ggalt::geom_encircle

#All time-points
pseq <- Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

p_legend <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) + # Matching colors for ellipses and points
  ggtitle("All time-points") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
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
pseq <- Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% c("Day 01"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p1 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) + # Matching colors for ellipses and points
  ggtitle("1 dpf") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
p1

#3 dpf
pseq <- Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% c("Day 03"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p3 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment),  shape = 24, size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) + 
  ggtitle("3 dpf") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
p3

#6 dpf
pseq <- Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% c("Day 06"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p6 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), shape = 22, size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) + 
  ggtitle("6 dpf") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
p6

#15 dpf
pseq <- Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% c("Day 15"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p15 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), shape = 3, size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#C77CFF", "lightgreen")) + 
  ggtitle("15 dpf") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
p15

#Spat
pseq <- Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% c("Spat"))
pseq.rel <- microbiome::transform(pseq, "compositional")


#note: for MU42022 no Heat treatment by Spat (blue colour)

ord <- ordinate(pseq.rel, "MDS", "bray")

pSpat <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), shape = 7, size = 6) +
  scale_colour_manual(values = c("#F8766D", "#C77CFF", "lightgreen")) +
  scale_fill_manual(values = c("#F8766D", "#C77CFF", "lightgreen")) + 
  ggtitle("Spat") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
pSpat

#arrange PCO plots -----

grid.arrange(p1, p3, p6, p15, pSpat, p, ncol = 3)

ggarrange(p1, p3, p6, p15, pSpat, p, ncol = 6, common.legend = TRUE, legend="bottom")

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


#mb2021 project ----

#Trying different ordination methods and distance methods
#Above code uses MDS/PCOa plots with Brays-Curtis distance matrix (relative abundance)


#NMDS = Performs Non-metric MultiDimenstional Scaling of a sample-wise ecological distance matrix 
#onto a user-specified number of axes, k.

pseq <- Marissa_mb2021_filtered_20240203
pseq <- subset_samples(pseq, Age %in% c("Spat"))
pseq <- subset_samples(pseq, !Family %in% c("9"))

pseq <- subset_samples(pseq, !Family %in% c("1"))

#correct family column for mb2021

pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

pseq@sam_data$Family <- as.character(pseq@sam_data$Family)

pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")


p <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Family") +
  geom_point(aes(fill = Treatment), size = 5) +
  scale_colour_manual(values = c("#F8766D", "#00ab41", "#1B98E0FF")) +
  ggtitle("Spat") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
p

p1 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), shape = 16, size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00ab41", "#1B98E0FF")) +
  ggtitle("1 dpf") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
p1


p2 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), shape = 17, size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00ab41", "#1B98E0FF")) +
  ggtitle("18 dpf") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
p2

p3 <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), shape = 15, size = 6) +
  scale_colour_manual(values = c("#F8766D", "#00ab41", "#1B98E0FF")) + 
  ggtitle("Spat") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
p3


grid.arrange(p1, p2, p3, p, ncol = 2)

ggarrange(p, p1, p2, p3, nrow = 2, ncol =2, common.legend = TRUE, legend="bottom")

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




#want to try weighted and unweighted Unifrac when I have the phylo tree (need raw seqs to get)
