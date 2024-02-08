#### Inval Tutorial Draft ----

### Source:  https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html 


#to look for differentially abundant ASVs using indicator species analysis


#Packages ----

install.packages("indicspecies")
library(indicspecies)
library(ggplot2)
#Set data ----

Marissa_MU42022_rare <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rare.rds")

pseq <- Marissa_MU42022_rare

#Load objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree


#Extract abundance matrix ----
#from the phyloseq object using phyloseq

OTU1 = as(OTU, "matrix")
write.csv(OTU1, file="Data_fram_1.cvs",row.names=TRUE)

write.table(OTU1,file="data_table.csv",sep=",",dec = " ")
####Format to example data and reload below for actual test 

#reload edited table
data_table <- read.csv("data_table.csv")

pc_FUN = read.csv("data_table.csv", header= TRUE)

pc_FUN <- data_table

#if removing samples ----

#Day 1 only ----

pc_FUN <- data_table[data_table$Age == "Day 01", ]


#Spat only ----

pc_FUN <- data_table[data_table$Age == "Spat", ]


#Larvae only ----

pc_FUN <- data_table[data_table$Age %in% c("Day 01", "Day 03", "Day 06", "Day 15"), ]

  
####Test ASVs ----

#Inverse data
funi_df<- t(pc_FUN)

###make into a matrix and populate::: This tells r what is metadata and what is the actual data ... Below 5-952 are the coloumns that are the data

matrix_F = pc_FUN[ ,8:366]

### Make the equation. Saying we want to examine specific column of metadata
time_a_F = pc_FUN$Treatment

### Run test 
inv_F = multipatt(matrix_F, time_a_F, func = "r.g", control = how(nperm=9999))
results <- summary(inv_F)

###Example of using the data data but testing different variables (from metadata)

trt_a_F = pc_FUN$Treatment
Just_trt_inv_F = multipatt(matrix_F, trt_a_F, func = "r.g", control = how(nperm=9999))
summary(Just_trt_inv_F)

#subset data

#Bubble plot results ----
#ASV7 and 18 unique to PB/PBH treatments
#plotting ASV 7 and 18 on day 1 data as was differentially abundant on day 1 (p < 0.05)

pc_FUN <- day1_PB_ASV

pc_FUN <- pc_FUN[pc_FUN$Age == "Day 01", ]


ggplot(pc_FUN, aes(x=Treatment, y= ASV, size = Abundance, color = Treatment)) + 
  geom_point(alpha=0.7)+ 
  scale_size(range = c(.1, 10)) +
  scale_colour_ipsum() +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("") +
  xlab("") +
  theme(legend.position = "none")

#plotting ASV 376 and 88 on spat data as was differentially abundant in spat (p < 0.05)

pc_FUN <- PB_ASVs

pc_FUN <- pc_FUN[pc_FUN$Age == "Spat", ]

pc_FUN <- pc_FUN[pc_FUN$ASV %in% c("ASV88", "ASV376"), ]


ggplot(pc_FUN, aes(x=Treatment, y= Genus, size = Abundance, color = Treatment)) + 
  geom_point(alpha=0.7)+ 
  scale_size(range = c(.1, 10)) +
  scale_colour_ipsum() +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("") +
  xlab("") +
  theme(legend.position = "none")