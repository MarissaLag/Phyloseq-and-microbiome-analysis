#### Inval Tutorial Draft ----

### Source:  https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html 

#Packages ----

install.packages("indicspecies")
library(indicspecies)
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
data_table <- read_csv("data_table.csv")

pc_FUN = read.csv("data_table.csv", header= TRUE)


####Test ASVs ----

#Inverse data
funi_df<- t(pc_FUN)

###make into a matrix and populate::: This tells r what is metadata and what is the actual data ... Below 5-952 are the coloumns that are the data

matrix_F = pc_FUN[ ,8:366]

### Make the equation. Saying we want to examine specific column of metadata
time_a_F = pc_FUN$TreatmentxAge

### Run test 
inv_F = multipatt(matrix_F, time_a_F, func = "r.g", control = how(nperm=9999))
summary(inv_F)

###Example of using the data data but testing different variables (from metadata)

trt_a_F = pc_FUN$Treatment
Just_trt_inv_F = multipatt(matrix_F, trt_a_F, func = "r.g", control = how(nperm=9999))
summary(Just_trt_inv_F)
