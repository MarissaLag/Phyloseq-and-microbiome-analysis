##take out Spat samples from data - rerun PCOA-MDS script
#as well, later on, take out F2,3,4 and only analyze F1 data


##### 00. Front Matter #####
# Clear space
# rm(list=ls())

## Load libraries
library("devtools")
library(phyloseq)
library(microbiome)
library(ggalt)
library(ggplot2)

## setwd
setwd("C:/Users/maris/OneDrive/Documents/USRA2021/mb2021/Data")

##### 01. Load data #####
#data <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022.rds")
pseq <- Marissa_MU42022_rarefied_20231016
pseq<- Marissa_mb2021_filtered_20240203

# Different project
#Marissa_Oyster <- Rare_filtered_data

#theme

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

#remove samples
#Remove day 3 (only 1 sample remaining) for mb2021 project

pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))

#change data genetics into character

pseq@sam_data$Genetics <- as.character(pseq@sam_data$Genetics)

str(pseq@sam_data)

##### 02. Analysis #####
head(pseq@sam_data)
str(pseq@sam_data)
unique(pseq@sam_data$Age)

#ALL samples ----
#reset pseq

#convert to compositional data

pseq.rel <- microbiome::transform(pseq, "compositional")

pseq_fam <- microbiome::aggregate_rare(pseq,level = "Family", detection = 50/100, prevalence = 70/100)

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq.rel, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq, ord, color = "Treatment", shape = "Age") + geom_point(size = 4)

plot_ordination(pseq, ord, color = "Treatment", shape = "Age") +
  geom_point(size = 3) +
  scale_colour_brewer(palette="Set1") +
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 11, face = "plain"), 
        axis.title = element_text(size = 13)) +
  ggalt::geom_encircle(aes(fill = Treatment), color = "black", expand = 0.2, alpha = 0.2)
  

?ggalt::geom_encircle



#All larvae ----

# "tolower" fxn can be used to convert all Spat to spat (all uppercase to lowercase)

pseq
pseq.nospat <- subset_samples(pseq, !(Age %in% c("Spat")))
pseq.nospat

#also remove algae

pseq.nospat <- subset_samples(pseq.nospat, !(Sample.type %in% c("Algae")))

#composition
pseq_fam <- microbiome::aggregate_rare(pseq.nospat, level = "Family", detection = 50/100, prevalence = 70/100)

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq.core, "MDS", "bray")

plot_ordination(pseq.nospat, ord, color = "Treatment", shape = "Age") +
  geom_point(size = 4)

plot_ordination(pseq.core, ord, color = "Treatment", shape = "Age") +
  geom_point(size = 4) +
  ggalt::geom_encircle(aes(fill = Treatment), color = "black", expand = 0.2, alpha = 0.2)

p = p + geom_point(size=7, alpha=0.75)
p = p + scale_colour_brewer(type="qual", palette="Set1")
p

#Spat only ----
#pseq <- Marissa_MU42022_rarefied_20231016

pseq<- Marissa_mb2021_filtered_20240203

pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))


#Compositional ----

pseq.rel <- microbiome::transform(pseq, "compositional")

#pseq <- subset_samples(pseq, !(Age %in% c("Day 03", "Day 06", "Day 15", "Day 01")))

#pseq <- subset_samples(pseq, !(Sample.type %in% "Algae"))


set.seed(4235421)

ord <- ordinate(pseq.rel, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(size = 3) +
  scale_colour_brewer(palette="Set1") +
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 11, face = "plain"), 
        axis.title = element_text(size = 13)) +
  ggalt::geom_encircle(aes(fill = Treatment), color = "black", expand = 0.2, alpha = 0.2)




# plot MDS/PcoA

p <- plot_ordination(pseq.rel, ord, color = "Treatment", label = "Library_Name") + geom_point(size = 4)
p <- p + scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF"))
p <- p + ggtitle("18 dpf")
p <- p + theme(plot.title = element_text(hjust = 0.5))  
p4 <- p + ggalt::geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)  

print(p4)


#Make plot for each time-point

#pseq<- Marissa_mb2021_filtered_20240203
#pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))

pseq <- Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

pseq.rel <- microbiome::transform(pseq, "compositional")

pseq.rel <- subset_samples(pseq.rel, Age %in% "Day 01")

set.seed(4235421)

ord <- ordinate(pseq.rel, "MDS", "bray")

p <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") + geom_point(size = 4)
p <- p + scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF","lightgreen"))
p <- p + ggtitle("1 dpf")
p <- p + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="none")
p1 <- p + ggalt::geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)  

print(p1)

p <- plot_ordination(pseq.rel, ord, color = "Treatment") + geom_point(size = 4, shape = 17)
p <- p + scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF"))
p <- p + ggtitle("18 dpf")
p <- p + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="none")
p2 <- p + ggalt::geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)

print(p2)

p <- plot_ordination(pseq.rel, ord, color = "Treatment") + geom_point(size = 4, shape = 15 )
p <- p + scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF"))
p <- p + ggtitle("Spat")
p <- p + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="none")
p3 <- p + ggalt::geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2) 

print(p3)

p <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age", label = "Library_Name") + geom_point(size = 4)
p <- p + scale_colour_manual(values = c("#F8766D", "#00BFC4", "#C77CFF"))
p <- p + ggtitle("All Time-points")
p <- p + theme(plot.title = element_text(hjust = 0.5))
p4 <- p + ggalt::geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2) 

print(p4)


#arrange pco plots

library(gridExtra)
library(ggpubr)

grid.arrange(p1, p2, p3, p4, ncol = 2)

ggarrange(p4, p1, p2, p3, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p4)

p_combined <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               p3 + theme(legend.position="none"),
                               p4 + theme(legend.position="none"),
                               ncol=2),
                   mylegend, nrow=1)


## Answer key colour ----
# scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"))


#DAY 1 ONLY, no algae ----

#reload data

pseq <- Marissa_MU42022_rarefied_20231016

pseq <- subset_samples(pseq, !(Age %in% c("Day 03", "Day 06", "Day 15", "Spat")))

pseq <- subset_samples(pseq, !(Sample.type %in% "Algae"))


pseq_fam <- microbiome::aggregate_rare(pseq,level = "Family", detection = 50/100, prevalence = 70/100)
#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq.core, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq.core, ord, color = "Treatment") + geom_point(size = 4) + theme_classic()

plot_ordination(pseq.core, ord, color = "Treatment") +
  geom_point(size = 4) +
  ggalt::geom_encircle(aes(fill = Treatment), color = "black", expand = 0.2, alpha = 0.2) +
  theme_classic()


#DAY 3 ONLY, no algae ----

#reload data

pseq <- Marissa_MU42022_rarefied_20231016

pseq <- subset_samples(pseq, !(Age %in% c("Day 01", "Day 06", "Day 15", "Spat")))

pseq <- subset_samples(pseq, !(Sample.type %in% "Algae"))


pseq_fam <- microbiome::aggregate_rare(pseq,level = "Family", detection = 50/100, prevalence = 70/100)
#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq.core, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq.core, ord, color = "Treatment") + geom_point(size = 4) + theme_classic()

plot_ordination(pseq.core, ord, color = "Treatment") +
  geom_point(size = 4) +
  ggalt::geom_encircle(aes(fill = Treatment), color = "black", expand = 0.2, alpha = 0.2) +
  theme_classic()



#DAY 6 ONLY, no algae ----

#reload data

pseq <- Marissa_MU42022_rarefied_20231016

pseq <- subset_samples(pseq, !(Age %in% c("Day 01", "Day 03", "Day 15", "Spat")))

pseq <- subset_samples(pseq, !(Sample.type %in% "Algae"))


pseq_fam <- microbiome::aggregate_rare(pseq,level = "Family", detection = 50/100, prevalence = 70/100)
#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq, ord, color = "Treatment") + geom_point(size = 4) + theme_classic()

plot_ordination(pseq, ord, color = "Treatment") +
  geom_point(size = 4) +
  ggalt::geom_encircle(aes(fill = Treatment), color = "black", expand = 0.2, alpha = 0.2) +
  theme_classic()

#DAY 15 ONLY, no algae ----

#reload data

pseq <- Marissa_MU42022_rarefied_20231016

pseq <- subset_samples(pseq, !(Age %in% c("Day 01", "Day 03", "Day 06", "Spat")))

pseq <- subset_samples(pseq, !(Sample.type %in% "Algae"))


pseq_fam <- microbiome::aggregate_rare(pseq,level = "Family", detection = 50/100, prevalence = 70/100)
#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq.core, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq.core, ord, color = "Treatment") + geom_point(size = 4) + theme_classic()

plot_ordination(pseq.core, ord, color = "Treatment") +
  geom_point(size = 4) +
  ggalt::geom_encircle(aes(fill = Treatment), color = "black", expand = 0.2, alpha = 0.2) +
  theme_classic()







#below code for Marissa MB2021 project

##DAY 18 SAMPLES ONLY ----

pseq_filtered <- subset_samples(Marissa_Oyster, !(Age %in% c("Spat", "Day 1", "Day 3", "Day 6")))

pseq_filtered <- subset_samples(Marissa_Oyster, !(Sample.type %in% c("Algae")))


pseq_fam <- microbiome::aggregate_rare(pseq_filtered, level = "Family", detection = 50/100, prevalence = 70/100)
#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq_filtered, "MDS", "bray")

plot_ordination(pseq_filtered, ord, color = "Tank_treatment", shape = "Age") +
  geom_point(size = 4) +
  ggalt::geom_encircle(aes(fill = Tank_treatment), color = "black", expand = 0.2, alpha = 0.2) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), labels = c("Control", "High salinity", "Low salinity")) +
  scale_shape_manual(values = c("day_18" = 17), labels = c("Day 18")) + theme_bw()

##remove DAY 3, day 1, and day 18 = sPAT SAMPLES ONLY ----

pseq_filtered <- subset_samples(Marissa_Oyster, !(Age %in% c("day_18", "day_1")))

pseq_fam <- microbiome::aggregate_rare(pseq_filtered, level = "Family", detection = 50/100, prevalence = 70/100)
#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq_filtered, "MDS", "bray")


plot_ordination(pseq_filtered, ord, color = "Tank_treatment", shape = "Age") +
  geom_point(size = 4) +
  ggalt::geom_encircle(aes(fill = Tank_treatment), color = "black", expand = 0.2, alpha = 0.2) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), labels = c("Control", "High salinity", "Low salinity")) +
  scale_shape_manual(values = c("spat" = 3), labels = c("Spat")) + theme_bw()

##all samples PCO with ellipses ----



#create objects

pseq <- Marissa_Oyster

pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq, "MDS", "bray")

#scree plot - to assess pca axes

eigenvalues <- ord$values$Eigenvalues
eigenvalues <- eigenvalues^2

#plot MDS/PcoA


plot_ordination(pseq, ord, color = "Tank_treatment", shape = "Age") +
  geom_point(size = 4) +
  ggalt::geom_encircle(aes(fill = Tank_treatment), color = "black", expand = 0.2, alpha = 0.2) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), labels = c("Control", "High salinity", "Low salinity")) +
  scale_shape_manual(values = c("day_1" = 16, "day_18" = 17, "day_3" = 15, "spat" = 3), 
                     labels = c("Day 1", "Day 18", "Day 3", "Spat")) + theme_bw()



plot_ordination(pseq, ord, color = "Tank_treatment", shape = "Age") +
  geom_point(size = 4) +
  ggalt::geom_encircle(aes(fill = Tank_treatment), color = "black", expand = 0.2, alpha = 0.2) +
  scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8"), 
                     labels = c("Control", "High salinity", "Low salinity")) +
  scale_shape_manual(values = c("day_1" = 16, "day_18" = 17, "day_3" = 15, "spat" = 3), 
                     labels = c("Day 1", "Day 18", "Day 3", "Spat")) + theme_bw()







##issue with comparing multiple plot ordinations = age groups have varying number of samples = we can't know if plot ordination is an artificate of the microbiome data itself OR number of samples

##try to rarify sample data and rarify abundance data (as I don't think it has been rarified number of sequences (abundance was rarified to >0.05%)? But not sure...)

# Assuming you have a phyloseq object named 'pseq' containing abundance data
# and 'sample_data' with sample metadata

# Load or install the 'phyloseq' package if you haven't already

library(phyloseq)

pseq <- Marissa_Osyter
OTU = Marissa_Osyter@otu_table


##first fix Genetics Metadata for spat samples if not done yet

# Replace "9" with "1" in group_samples vector



# Step 1: Create a function to rarify a given data matrix to a specified count
rarify_data <- function(data_matrix, count) {
  total_count <- sum(data_matrix)
  if (total_count <= count) {
    # If the total count is less than or equal to the desired count, return the data as is
    return(data_matrix)
  } else {
    # Otherwise, rarify the data to the desired count using the "rarefy" function from the vegan package
    library(vegan)
    rarified_data <- rarefy(data_matrix, sample = count)
    return(rarified_data)
  }
}

# Step 2: Group your data by the "Age" variable
library(dplyr)
grouped_data <- abundance_data %>%
  left_join(group_samples, by = "Library_Name") %>%
  group_by(Age)

# Step 3: Within each group, rarify the samples to 3 if the number of samples is greater than 3
grouped_data <- grouped_data %>%
  mutate(rarified_abundance = rarify_data(abundance_data, 3))

# Step 4: After rarifying the samples within each group, rarify the entire data matrix to 5000
total_rarified_abundance <- rarify_data(sum(grouped_data$rarified_abundance), 5000)

# Now you have the rarefied data matrix with each group rarified to 3 samples and the entire data rarified to 5000.







**********************************************************************************************************************************
  #removing F2,3,4 from analysis since F1 appears to be different (especially at spat stage) from other samples = may be blurring treatment effects
  
  
  # Remove F2,3,4 samples from the phyloseq object - watch, the subset samples fxn is case sensitive
  #remember - Metadata file Family column is wrong for spat samples - tank # listed rather than family
  # F1 = 1, 9, 13
  # F2 = 2, 10, 14
  # F3 = 3, 11, 15
  # F4 = 4, 12, 16
  
pseq_filtered <- subset_samples(Marissa_Osyter, !(Family %in% c("2", "10", "14", "3", "11", "15", "4", "12", "16")))



pseq_fam <- microbiome::aggregate_rare(pseq_filtered, level = "Family", detection = 50/100, prevalence = 70/100)
#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq_filtered, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq_filtered, ord, color = "Tank_treatment", shape = "Age") + geom_point(size = 4)

##remove F1

pseq_filtered <- subset_samples(Marissa_Osyter, !(Family %in% c("1", "9", "13")))
pseq_fam <- microbiome::aggregate_rare(pseq_filtered, level = "Family", detection = 50/100, prevalence = 70/100)
#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq_filtered, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq_filtered, ord, color = "Tank_treatment", shape = "Age") + geom_point(size = 4)


#remove F1 and larval samples

##remove F1

pseq_filtered <- subset_samples(Marissa_Osyter, !(Family %in% c("1", "9", "13")))

pseq_filtered <-  subset_samples(pseq_filtered, !(Age %in% c("day_1", "day_3", "day_18")))

##remove T10r3 - seems like outlier

pseq_filtered <-  subset_samples(pseq_filtered, !(Library_Name %in% "T10r3"))

pseq_fam <- microbiome::aggregate_rare(pseq_filtered, level = "Family", detection = 50/100, prevalence = 70/100)
#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq_filtered, "MDS", "bray")

#plot MDS/PcoA

plot_ordination(pseq_filtered, ord, color = "Tank_treatment", shape = "Age") + geom_point(size = 4) + geom_text(aes(label = Library_Name), vjust = 0, nudge_y = 0.1)
                                                                                                                




