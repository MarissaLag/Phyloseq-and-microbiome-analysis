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

#set theme.marissa as plot theme

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


#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#PERMANOVA ----
#source: https://microbiome.github.io/tutorials/PERMANOVA.html

#set seed

set.seed(452)

#if any columns have missing values (NA), must remove

pseq_filtered <- subset_samples(pseq, !Treatment %in% NA)
pseq_filtered <- subset_samples(pseq_filtered, !Genetics %in% NA)

#different method to remove NA
Meta <- subset(Metadata, !is.na(Treatment))
Meta <- meta(pseq.fam.rel)

#use same data used in PCO plot
#data must be in data frame
#convert to data frame

Data <- data.frame(sample_data(pseq_filtered))

#create matrix - 2 methods below

Bray_dist <- ordinate(pseq_filtered, method = "MDS", distance = "bray", weighted = TRUE)

Bray_dist<- phyloseq::distance(pseq_filtered, method = "bray", weighted = TRUE)


#Run PERMANOVA
#data must be data frame (independent variables, i.e., your metadata)

permanova <- adonis2(Bray_dist ~ Treatment*Age*Genetics,
                     data = Data, permutations=999, method = "bray")
permanova

#For some reason, get different p values with only Treatment?

permanova2 <- adonis2(Bray_dist ~ Treatment,
                     data = Data, permutations=999, method = "bray")

permanova2

#below is for "strata" analysis - if factors within data are nested use strata
#I am still unsure when strata is needed - most examples are for repeated measures

permanova_strata <- adonis2(Bray_dist ~ Genetics,
                     data = Data, permutations=999, method = "bray", strata = pseq_filtered@sam_data$Treatment)

#make new data frame where treatment and age are numeric

Data$Treatment <- sub("Control", 1, Data$Treatment)
Data$Treatment <- sub("Probiotics.+.HT", 3, Data$Treatment)
Data$Treatment <- sub("Probiotics", 2, Data$Treatment)
Data$Treatment <- sub("High temperature", 4, Data$Treatment)

Data$Age <- sub("Day 01", 1, Data$Age)
Data$Age <- sub("Day 03", 3, Data$Age)
Data$Age <- sub("Day 06", 6, Data$Age)
Data$Age <- sub("Day 15", 15, Data$Age)
Data$Age <- sub("Spat", 20, Data$Age)

View(Data)

#making new data frame - code not working

replications <- c(11, 13, 11, 12, 20)

dat <- expand.grid(Age=gl(5, replications), Treatment=factor(c(1,2,3,4)),Genetics=factor(c(1,2,3,4)) )
View(dat)

#PCO plot ----

set.seed(4235421)

ord <- ordinate(pseq, "MDS", "bray")

#plot MDS/PcoA ----

plot_ordination(pseq, ord, color = "Treatment", shape = "Age") + geom_point(size = 4) + theme.marissa()

#to subset samples (e.g., time series analysis) and create more PCO plots see Factors_PCO_PCA_mb2021 tutorial

#scree plot to assess PCo axes see plot_ordination_methods script